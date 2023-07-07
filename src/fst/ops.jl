# SPDX-License-Identifier: CECILL-2.1

# Aliases to specify different behavior for acceptors and transducers.
const TransducerOrAcceptor = AbstractFST
const Transducer = AbstractFST{<:Semiring,<:LabelSymbolMapping}
const Acceptor = AbstractFST{<:Semiring,<:LabelSymbol}

function accessible(fst::AbstractFST)
	m = ones(Bool, nstates(fst))
	prune(xᵀ) = begin
		x = parent(xᵀ)
		vₙᵀ = transpose(x .* m)
		m[findnz(x)[1]] .= false
		vₙᵀ
	end

	T = S(fst) * transpose(D(fst))
	transpose(α(fst)) * powerseries(prune, T)
	.! m
end

function coaccessible(fst::AbstractFST)
	# Prune the state already visited
	m = ones(Bool, nstates(fst))
	prune(x) = begin
		x = x
		vₙ = x .* m
		m[findnz(x)[1]] .= false
		vₙ
	end

	T = S(fst) * transpose(D(fst))
	powerseries(prune, T) * ω(fst)
	.! m
end

#=====================================================================#
# Unary operations
#=====================================================================#

"""
    filterarcs(fn, fst)

Create a new FST by removing arcs `a = (s, d, l, w)` for which `fn(a)`
returns false.
"""
function filterarcs(fn::Function, fst::TransducerOrAcceptor)
	ii = findall(fn, collect(arcs(fst)))
	P = length(ii)
    M = sparse(ii, 1:P, one(semiring(fst)), narcs(fst), P)

	FST(α(fst), S(fst) * M, D(fst) * M, ω(fst), λ(fst)[ii])
end

"""
    filterstates(fn, fst)

Create a new FST by removing states `s = (q, iw, fw)` for which `fn(q)`
returns false. Arcs reaching and leaving `q` are also removed.
"""
function filterstates(fn::Function, fst::TransducerOrAcceptor)
	qq = findall(fn, collect(states(fst)))
	P = length(qq)
    Mᵀ = transpose(sparse(qq, 1:P, one(semiring(fst)), nstates(fst), P))

	f = filterarcs(fst) do (s, d, l, w)
		s ∈ Set(qq) && d ∈ Set(qq)
	end

    FST(Mᵀ * α(f), Mᵀ * S(f), Mᵀ * D(f), Mᵀ * ω(f), λ(f))
end

"""
    relabel(fn, fst)

Create a new FST where each label `l` is mapped `fn(l)`
"""
relabel(fn, fst::TransducerOrAcceptor) =
    FST(α(fst), S(fst), D(fst), ω(fst), map(fn, λ(fst)))

"""
    reverse(fst)

Create a new FST such that which accept the mirror image of every
string accepted by `fst`.
"""
Base.reverse(fst::TransducerOrAcceptor) =
    FST(ω(fst), D(fst), S(fst), α(fst), λ(fst))

"""
    inv(fst)

Create a new transducer where the transduction rules of `fst` (x => y)
are swapped (y => x).
"""
Base.inv(fst::Transducer) = relabel(l -> last(l) => first(l), fst)
Base.inv(fst::Acceptor) = fst

"""
    project(fst[; project_output = false])

Create an acceptor by copying every input label to its output label
or vice versa when `project_output = true`.
"""
project(fst::Transducer; project_output = false) =
    relabel(l -> project_output ? last(l) : first(l), fst)
project(fst::Acceptor; project_output = false) = fst

function connect(fst::TransducerOrAcceptor)
	qq = findall(identity, accessible(fst) .* coaccessible(fst))
	P = length(qq)
    M = sparse(qq, 1:P, one(semiring(fst)), nstates(fst), P)

	fst_pruned = filterarcs(fst) do (s, d, l, w)
		s ∈ Set(qq) && d ∈ Set(qq)
	end

	FST(
		transpose(M) * α(fst_pruned),
		transpose(M) * S(fst_pruned),
		transpose(M) * D(fst_pruned),
		transpose(M) * ω(fst_pruned),
		λ(fst_pruned)
	)
end

#=====================================================================#
# Binary operations
#=====================================================================#

function Base.kron(X::TransducerOrAcceptor, Y::TransducerOrAcceptor)
	FST(
		kron(α(X), α(Y)),
		kron(S(X), S(Y)),
		kron(D(X), D(Y)),
		kron(ω(X), ω(Y)),
		[(first(l1), first(l2)) => (last(l1), last(l2))
		 for l1 in λ(X) for l2 in λ(Y)]
	)
end

function _compose_ϵfree(match, A::AbstractFST, B::AbstractFST)
	# Tensor product between A and B
	AB = kron(A, B)

	# Remove arcs (lA, lB) for which the output label
	# of lA does not match the input label lB
	pruned_AB = filterarcs(match, AB)

	# Remove unreachable nodes the resulting FST
	result = connect(pruned_AB)

	# Transodrm the label (lA, lB)
	# into new labels lA[1]:lB[2]
	relabel(result) do l
		lA, lB = l
        i1, o1 = lA
        i2, o2 = lB
        o2 = o2 == -3  ? i2 : o2
		i1 => o2
	end
end

function _compose_ϵfree(A::Transducer, B::Transducer)
    _compose_ϵfree(A, B) do (s, d, l, w)
        lA, lB = l
        last(lA) == first(lB)
    end
end

function _compose(X::Transducer, Y::Transducer)
	ϵ, ϵ1, ϵ2, ϕ = 0, -1, -2, -3
	K = promote_type(semiring(X), semiring(Y))

	X = relabel(X) do l
		FiniteStateAutomata.isepsilon(last(l)) ? first(l) => ϵ2 : l
	end
	X̃ = FST(
		α(X),
		hcat(S(X), sparse(1:nstates(X), 1:nstates(X), one(semiring(X)))),
		hcat(D(X), sparse(1:nstates(X), 1:nstates(X), one(semiring(X)))),
		ω(X),
		vcat(λ(X), repeat([ϵ => ϵ1], nstates(X)))
	)

	Y = relabel(Y) do l
		FiniteStateAutomata.isepsilon(first(l)) ? ϵ1 => last(l) : l
	end
	Ỹ = FST(
		α(Y),
		hcat(S(Y), sparse(1:nstates(Y), 1:nstates(Y), one(semiring(Y)))),
		hcat(D(Y), sparse(1:nstates(Y), 1:nstates(Y), one(semiring(Y)))),
		ω(Y),
		vcat(λ(Y), repeat([ϵ2 => ϵ], nstates(Y)))
	)

	F = FST(
		sparsevec([1], one(K), 3),
		sparse([1, 1, 1, 1, 2, 2, 3, 3], 1:8, one(K), 3, 8),
		sparse([1, 1, 1, 1, 2, 2, 3, 3], [1, 2, 5, 7, 3, 6, 4, 8], one(K), 3, 8),
		sparsevec([1, 2, 3], one(K), 3),
		[ϕ => ϕ, ϵ2 => ϵ1, ϵ1 => ϵ1, ϵ2 => ϵ2, ϕ => ϕ, ϵ1 => ϵ1, ϕ => ϕ, ϵ2 => ϵ2]
	)

	XF = _compose_ϵfree(X̃, F) do (s, d, l, w)
        o1, i2 = first(last(l)), last(first(l))
        o1 == i2 || (i2 == ϕ && o1 > 0 )
	end

	_compose_ϵfree(XF, Ỹ) do (s, d, l, w)
        o1, i2 = first(last(l)), last(first(l))
        o1 == i2 || (o1 == ϕ && i2 > 0)
	end
end

function Base.:∘(X::Transducer, Y::Transducer)
    if any(l -> isoutputepsilon(l), λ(X)) || any(l -> isinputepsilon(l), λ(Y))
        _compose(X, Y)
    else
        _compose_ϵfree(X, Y)
    end
end

#"""
#    W(A)
#
#Return the total weight of `A`, i.e. the ``\\oplus``-sum of all the
#path's weight in `A`.
#"""
#W(A::AbstractFST) =
#    ρ(A) ⊕ (transpose(α(A)) * MatrixPowerSum(T(A)) * ω(A)) # TODO optimize with dot(., ., .)
#
#function closure(A::AbstractFST; plus = false)
#    K = eltype(α(A))
#    TA = T(A)
#    TB = TransitionMatrix(
#        TA.S,
#        hcat(TA.U, ω(A)),
#        blockdiag(TA.E, spzeros(K, 1, 1)),
#        vcat(TA.V, transpose(α(A)))
#    )
#    ρB = iszero(ρ(A)) && ! plus ? one(K) : ρ(A)
#    FST(α(A), TB, ω(A), ρB, λ(A))
#end
#
#"""
#    cat(A1[, A2, ...])
#
#Return the concatenation of the given FSTs.
#"""
#function Base.cat(A::AbstractFST, B::AbstractFST)
#    K = promote_type(eltype(α(A)), eltype(α(B)))
#    TA, TB = T(A), T(B)
#    TC = TransitionMatrix(
#        blockdiag(TA.S, TB.S),
#        hcat(blockdiag(TA.U, TB.U), vcat(ω(A), spzeros(K, nstates(B)))),
#        blockdiag(blockdiag(TA.E, TB.E), spzeros(K, 1, 1)),
#        vcat(blockdiag(TA.V, TB.V), transpose(vcat(spzeros(K, nstates(B)), α(B)))),
#    )
#
#    FST(
#        vcat(α(A), ρ(A) * α(B)),
#        TC,
#        vcat(ω(A) * ρ(B), ω(B)),
#        ρ(A) ⊗ ρ(B),
#        vcat(λ(A), λ(B))
#    )
#end
#Base.cat(A1::AbstractFST{K}, AN::AbstractFST{K}...) where K = foldl(cat, AN, init = A1)
#
#Π₁(A::Acceptor) = FST(α(A), T(A), ω(A), ρ(A), λ(A))
#Π₂(A::Acceptor) = FST(α(A), T(A), ω(A), ρ(A), λ(A))
#Π₁(A::AbstractFST) = FST(α(A), T(A), ω(A), ρ(A), first.(λ(A)))
#Π₂(A::AbstractFST) = FST(α(A), T(A), ω(A), ρ(A), last.(λ(A)))
#
#"""
#    union(A1[, A2, ...])
#    A1 ∪ A2
#
#Return the union of the given FST.
#"""
#function Base.union(A::AbstractFST, B::AbstractFST)
#    K = promote_type(eltype(α(A)), eltype(α(B)))
#    TA, TB = T(A), T(B)
#    TC = TransitionMatrix(
#        blockdiag(TA.S, TB.S),
#        blockdiag(TA.U, TB.U),
#        blockdiag(TA.E, TB.E),
#        blockdiag(TA.V, TB.V)
#    )
#
#    FST(vcat(α(A), α(B)), TC, vcat(ω(A), ω(B)), ρ(A) ⊕ ρ(B), vcat(λ(A), λ(B)))
#end
#
#Base.union(A1::AbstractFST{K}, AN::AbstractFST{K}...) where K = foldl(union, AN, init = A1)
