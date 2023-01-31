# SPDX-License-Identifier: CECILL-2.1

"""
    addskipedges(M, states, weights)

Added edges such that each state `i` in `states` is potentially
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

function _filter(A::AbstractFSA{K}, x::AbstractVector{Bool}) where K
    J = findall(x)
    M = sparse(1:length(J), J, one(K), length(J), nstates(A))
    FSA(M * α(A), M * T(A) * M', M * ω(A), ρ(A), λ(A)[J])
end

function Base.filter(f::Function, A::AbstractFSA)
    _filter(A, f.(1:nstates(A)))
end

"""
    notaccessible(A::AbstractFSA{K}) where K

Returns a vector `x` of dimension `nstates(A)` where `x[i]` is `one(K)`
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
coaccessible(A::AbstractFSA) = accessible(A |> reverse)

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
    FSA(
        Zₐ * α(A),
        Z .* T(A),
        Z[:, 1] .* ω(A),
        Zₐ * ρ(A),
        λ(A)
    )
end

function Base.replace(M::AbstractFSA, Ms)
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
        vcat([λ(M)[i] * λ(m) for (i, m) in enumerate(Ms)]...)
    )
end

function Base.replace(new::Function, M::AbstractFSA)
    replace(M, [new(i) for i in 1:nstates(M)])
end

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

"""
    globalrenorm(A::AbstractAcyclicFSA)

Global renormalization of the weights of `A`.
"""
function globalrenorm(A::AbstractAcyclicFSA)
    # Accumulate the weight backward starting from the end state.
    v = ω(A)
    σ = v
    while nnz(v) > 0
        v = T(A) * v
        σ += v
    end
    σ

    D = spdiagm(σ)
    FSA(σ .* α(A), T(A) * D, ω(A), ρ(A), λ(A)) |> renorm
end

