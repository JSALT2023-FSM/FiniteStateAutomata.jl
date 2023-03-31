# SPDX-License-Identifier: CECILL-2.1

"""
    addskipedges(M, states, weights)

Add edges such that each state `i` in `states` is potentially
skipped with weight `weights[i]`. This operation does not remove
any existing edge.
"""
function addskipedges(M, states, weights)
	Tₙ = T(M)
	αₙ = α(M)
	ωₙ = ω(M)
	ρₙ = ρ(M)
	Q = nstates(M)

	for (i, w) in zip(states, weights)
		σeeᵀ = sparse([i], [i], [w], Q, Q)
		αₙ = αₙ + Tₙ' * σeeᵀ * αₙ
		Tₙ = Tₙ + Tₙ * σeeᵀ * Tₙ
		ωₙ = ωₙ + Tₙ * σeeᵀ * ωₙ
		ρₙ = ρₙ + αₙ' * σeeᵀ * ωₙ
	end

	FSA(αₙ, Tₙ, ωₙ, ρₙ, λ(M))
end

"""
    Base.sum(A[; init = initstates(A), n = nstates(A)])

Accumulate the weight of all the paths of length less or equal to
`n`.
"""
function Base.sum(A::AbstractFSA; init = α(A), n = nstates(A) + 2)
    v = init
    σ = dot(v, ω(A)) + ρ(A)
    for i in 2:n-1
        v = T(A)' * v
        σ += dot(v, ω(A))
    end
    σ
end


"""
    cat(A1[, A2, ...])

Return the concatenation of the given FSA.
"""
function Base.cat(A1::AbstractFSA{K}, A2::AbstractFSA{K}) where K
    D1, D2 = size(T(A1), 1), size(T(A2), 1)
    FSA(
        vcat(α(A1), iszero(ρ(A1)) ? spzeros(K, D2) :  ρ(A1) * α(A2)),
        [T(A1) (ω(A1) * α(A2)');
         spzeros(K, D2, D1) T(A2)],
        vcat(iszero(ρ(A2)) ? spzeros(K, D1) : ω(A1) * ρ(A2), ω(A2)),
        ρ(A1) * ρ(A2),
        vcat(λ(A1), λ(A2))
    )
end
Base.cat(A1::AbstractFSA{K}, AN::AbstractFSA{K}...) where K =
    foldl(cat, AN, init = A1)

"""
    closure(A; plus = false)

Return the closure (or the closure plus if `plus` is `true`) of the FSA.
"""
function closure(A::AbstractFSA{K}; plus = false) where K
    FSA(
        α(A),
        T(A) + ω(A) * α(A)',
        ω(A),
        iszero(ρ(A)) && ! plus ? one(K) : ρ(A),
        λ(A)
    )
end

function _fsadet(M, edges, states)
	state2idx = Dict(q => i for (i, q) in enumerate(sort(collect(states))))
    D, Q = length(states), length(λ(M))
    K = eltype(α(M))

	I, V = [], K[]
	for i in edges[()]
		eᵢ = sparsevec(collect(i), one(K), Q)
		push!(I, state2idx[i])
        push!(V, dot(α(M), eᵢ))
	end
	a = sparsevec(I, V, D)

	I, J, V = [], [], K[]
	for i in keys(edges)
		i == () && continue
		for j in edges[i]
			push!(I, state2idx[i])
			push!(J, state2idx[j])
			eᵢ = sparsevec(collect(i), one(K), Q)
			eⱼ = sparsevec(collect(j), one(K), Q)
            push!(V, eᵢ' * T(M) * eⱼ)
		end
	end
	B = sparse(I, J, V, D, D)

	I, V = [], K[]
	for i in states
		eᵢ = sparsevec(collect(i), one(K), Q)
        v = dot(ω(M), eᵢ)
        if v ≠ zero(K)
            push!(I, state2idx[i])
            push!(V, v)
        end
	end
	o = sparsevec(I, V, D)

    l = [λ(M)[q[1]] for q in sort(collect(states))]

    FSA(a, B, o, ρ(M), l)
end

function _spmap(K, mapping)
	unique_ids = Set(mapping)
	id2idx = Dict(id => i for (i, id) in enumerate(collect(unique_ids)))
	idx_mapping = [id2idx[id] for id in mapping]

	Q, D = length(mapping), length(unique_ids)
	sparse(1:Q, idx_mapping, one(K), Q, D)
end

"""
    determinize(A::AbstractFSA)

Return an "equivalent" deterministic version of `A`. Note that the
weights of the merged paths are simply added and locally renormalized.
Therefore, the resulting FSA may have different weighting than the
original FSA.
"""
function determinize(A::AbstractFSA{K}) where K
    C = _spmap(K, λ(A))

	# Edges of the determinized FSA.
	edges = Dict()

	# Visited states of the determinized FSA.
	visited = Set()

	# Initialize the stack, `()` represents the initial state
	# of the determinized FSA.
    stack = Any[((), α(A))]
	while ! isempty(stack)
		# Get the first element from the stack.
		# `a` is the ancestor state (of the new FSA)
		# and `v` is the children of `a` in `A`.
		a, v = pop!(stack)

		# Merge states with the same label.
		V = spdiagm(v) * C

		# Iterate over the column of V that have non-zero elements.
		for i in 1:size(V, 2)
			c = V[:, i]
			nnz(c) == 0 && continue

			# A new state is defined based on the indices
			# of the non-zero elements.
			q = Tuple(findnz(c)[1])

			# Record the edge between the ancestor and the new
			# state.
			edges[a] = push!(get(edges, a, []), q)

			# If the new state `q` has not been visited yet
			# add its children (represented as a vector)
			# on the stack.
			if q ∉ visited
                push!(stack, (q, T(A)' * c))
				push!(visited, q)
			end
		end
	end
    _fsadet(A, edges, visited)
end

"""
    minimize(A::AbstractFSA)

Return a minimal equivalent FSA.
"""
minimize(A::AbstractFSA) = (reverse ∘ determinize ∘ reverse ∘ determinize)(A)

function _filter(A::AbstractFSA{K}, x::AbstractVector{Bool}) where K
    J = findall(x)
    M = sparse(1:length(J), J, one(K), length(J), nstates(A))
    FSA(M * α(A), M * T(A) * M', M * ω(A), ρ(A), λ(A)[J])
end

function Base.filter(f::Function, A::AbstractFSA)
    _filter(A, f.(1:nstates(A)))
end

"""
    accessible(A::AbstractFSA{K}) where K

Return a vector `x` of dimension `nstates(A)` where `x[i]` is `one(K)`
if the state `i` is not accessible, `zero(K)` otherwise.
"""
function accessible(A::AbstractFSA{K}) where K
    vₙ = α(A)

    m = ones(Bool, nstates(A))
	m[findnz(vₙ)[1]] .= false

	while nnz(vₙ) > 0
		uₙ = T(A)' * vₙ
		vₙ = uₙ .* m
		m[findnz(uₙ)[1]] .= false
	end
	.! m
end

"""
    coaccessible(A::AbstractFSA{K}) where K

Returns a vector `x` of dimension `nstates(A)` where `x[i]` is `one(K)`
if the state `i` is not accessible, `zero(K)` otherwise.
"""
function coaccessible(A::AbstractFSA)
    vₙ = ω(A)

    m = ones(Bool, nstates(A))
	m[findnz(vₙ)[1]] .= false

	while nnz(vₙ) > 0
		uₙ = T(A) * vₙ
		vₙ = uₙ .* m
		m[findnz(uₙ)[1]] .= false
	end
	.! m
end

"""
    connect(A::AbstractFSA)

Return a FSA ``B`` equivalent of ``A`` such that all states in ``B``
are accessible and coaccessible.
"""
function connect(A::AbstractFSA{K}) where K
    m = accessible(A) .* coaccessible(A)
    I = findall(m)
    M = sparse(I, 1:length(I), one(K), nstates(A), length(I))
    FSA(
        M' * α(A),
        M' * T(A) * M,
        M' * ω(A),
        ρ(A),
        λ(A)[I]
    )
end

"""
    renorm(A::FSA)

Local renormalization. For a state q, multiply each arc's weight by the
inverse of the sum of all the arcs' weight leaving q. In the resulting
FSA, sum of the all the arcs'weight of a state sum up to the semiring
one.
"""
function renorm(A::AbstractFSA)
    Z = inv.((sum(T(A), dims=2) .+ ω(A)))
    Zₐ = inv(sum(α(A)) + ρ(A))
    typeof(A)(
        Zₐ * α(A),
        Z .* T(A),
        Z[:, 1] .* ω(A),
        Zₐ * ρ(A),
        λ(A)
    )
end

@inline mergelabels(x) = x
@inline mergelabels(x, y) = (x..., y...)
@inline mergelabels(x, y, z...) = mergelabels(mergelabels(x, y), z...)

function Base.replace(f::Function, M::AbstractFSA, Ms)
    R = ρ.(Ms)
    I = findall(! iszero, R)
    M = addskipedges(M, I, R[I])

    A = blockdiag([α(m)[:,1:1] for m in Ms]...)
    Ω = blockdiag([ω(m)[:,1:1] for m in Ms]...)
    D = spdiagm([ρ(m) for m in Ms])
    FSA(
        vcat([α(M)[i] * α(m) for (i, m) in enumerate(Ms)]...),
        blockdiag([T(m) for m in Ms]...) + Ω * (T(M) + T(M) * D * T(M)') * A',
        vcat([ω(M)[i] * ω(m) for (i, m) in enumerate(Ms)]...),
        ρ(M),
        vcat([f.([λ(M)[i]], λ(m)) for (i, m) in enumerate(Ms)]...)
    )
end

function Base.replace(new::Function, M::AbstractFSA, labelfn::Function)
    replace(labelfn, M, [new(i) for i in 1:nstates(M)])
end

Base.replace(new::Function, M::AbstractFSA) =
    replace(new, M, mergelabels)

"""
    reverse(A)

Return the reversal of A.
"""
Base.reverse(A::AbstractFSA) = typeof(A)(ω(A), copy(T(A)'), α(A), ρ(A), λ(A))

"""
    union(A1[, A2, ...])
    A1 ∪ A2

Return the union of the given FSA.
"""
function Base.union(A1::AbstractFSA{K}, A2::AbstractFSA{K}) where K
    FSA(
        vcat(α(A1), α(A2)),
        blockdiag(T(A1), T(A2)),
        vcat(ω(A1), ω(A2)),
        ρ(A1) + ρ(A2),
        vcat(λ(A1), λ(A2))
    )
end
Base.union(A1::AbstractFSA{K}, AN::AbstractFSA{K}...) where K =
    foldl(union, AN, init = A1)

#= Functions for acyclic FSA =#

Base.cat(A1::AbstractAcyclicFSA, A2::AbstractAcyclicFSA) =
    AcyclicFSA(cat(parent(A1), parent(A2)))
minimize(A::AbstractAcyclicFSA) = AcyclicFSA(parent(A) |> minimize)
renorm(A::AbstractAcyclicFSA) = AcyclicFSA(parent(A) |> renorm)
Base.reverse(A::AbstractAcyclicFSA) = AcyclicFSA(parent(A) |> reverse)
Base.union(A1::AbstractAcyclicFSA, A2::AbstractAcyclicFSA) =
    AcyclicFSA(union(parent(A1), parent(A2)))

"""
    propagate(A::AbstractAcyclicFSA[; forward = true])

Multiply the weight of each arc by the sum-product of all the prefix
paths. If `forward` is `false` the arc's weight is multiply by the
sum-product of all the suffix paths of this arc.
"""
function propagate(A::AbstractAcyclicFSA; forward = true)
    v = forward ? α(A) : ω(A)
    M = forward ? T(A)' : T(A)
    σ = v
    while nnz(v) > 0
        v = M * v
        σ += v
    end
    σ

    D = spdiagm(σ)
    AcyclicFSA(FSA(σ .* α(A), T(A) * D, ω(A), ρ(A), λ(A)))
end

"""
    globalrenorm(A::AbstractAcyclicFSA)

Global renormalization of the weights of `A`.
"""
globalrenorm(A::AbstractAcyclicFSA) = propagate(A; forward = false) |> renorm

