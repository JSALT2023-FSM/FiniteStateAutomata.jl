# SPDX-License-Identifier: CECILL-2.1

"""
    Base.cumsum(A[; init = initstates(A), n = nstates(A)])

Accumulate the weight of all the paths of length less or equal to
`n`.
"""
function Base.sum(A::AbstractFSA; init = α(A), n = nstates(A))
    v = init
    σ = dot(v, ω(A)) + ρ(A)
    for i in 1:n
        v = T(A)' * v
        σ += dot(v, ω(A))
    end
    σ
end


"""
    union(A1[, A2, ...])
    A1 ∪ A2

Return the union of the given FSA.
"""
function Base.union(A1::AbstractFSA{K,L}, A2::AbstractFSA{K,L}) where {K,L}
    FSA(
        vcat(α(A1), α(A2)),
        blockdiag(T(A1), T(A2)),
        vcat(ω(A1), ω(A2)),
        ρ(A1) + ρ(A2),
        vcat(λ(A1), λ(A2))
    )
end
Base.union(A1::AbstractFSA{K,L}, AN::AbstractFSA{K,L}...) where {K,L} =
    foldl(union, AN, init = A1)

"""
    cat(A1[, A2, ...])

Return the concatenation of the given FSA.
"""
function Base.cat(A1::AbstractFSA{K,L}, A2::AbstractFSA{K,L}) where {K,L}
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
Base.cat(A1::AbstractFSA{K,L}, AN::AbstractFSA{K,L}...) where {K,L} =
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

"""
    reverse(A)

Return the reversal of A.
"""
Base.reverse(A::AbstractFSA) = FSA(ω(A), copy(T(A)'), α(A), ρ(A), λ(A))

"""
    renorm(A::FSA)

Local renormalization. For a state q, multiply each arc's weight by the
inverse of the sum of all the arcs' weight leaving q. In the resulting
FSA, sum of the all the arcs'weight of a state sum up to the semiring
one.
"""
function renorm(A::FSA{K}) where K
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

