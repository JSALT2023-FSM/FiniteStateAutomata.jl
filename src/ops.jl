# SPDX-License-Identifier: CECILL-2.1

"""
    union(fsa1[, fsa2, ...])
    fsa1 ∪ fsa2

Return the union of the given FSA.
"""
function Base.union(fsa1::AbstractFSA{K,L}, fsa2::AbstractFSA{K,L}) where {K,L}
    FSA(
        vcat(α(fsa1), α(fsa2)),
        blockdiag(T(fsa1), T(fsa2)),
        vcat(ω(fsa1), ω(fsa2)),
        ρ(fsa1) + ρ(fsa2),
        vcat(λ(fsa2), λ(fsa2))
    )
end
Base.union(fsa1::AbstractFSA{K,L}, fsaN::AbstractFSA{K,L}...) where {K,L} =
    foldl(union, fsaN, init = fsa1)

"""
    cat(fsa1[, fsa2, ...])

Return the concatenation of the given FSA.
"""
function Base.cat(fsa1::AbstractFSA{K,L}, fsa2::AbstractFSA{K,L}) where {K,L}
    D1, D2 = size(T(fsa1), 1), size(T(fsa2), 1)
    FSA(
        vcat(α(fsa1), zero(K) * α(fsa2)),
        [T(fsa1) (ω(fsa1) * α(fsa2)');
         spzeros(K, D2, D1) T(fsa2)],
        vcat(zero(K) * ω(fsa1), ω(fsa2)),
        ρ(fsa1) * ρ(fsa2),
        vcat(λ(fsa2), λ(fsa2))
    )
end
Base.cat(fsa1::AbstractFSA{K,L}, fsaN::AbstractFSA{K,L}...) where {K,L} =
    foldl(cat, fsaN, init = fsa1)

"""
    closure(fsa; plus = false)

Return the closure (or the closure plus if `plus` is `true`) of the FSA.
"""
function closure(fsa::AbstractFSA{K}; plus = false) where K
    FSA(
        α(fsa),
        T(fsa) + ω(fsa) * α(fsa)',
        ω(fsa),
        iszero(ρ(fsa)) && ! plus ? one(K) : ρ(fsa),
        λ(fsa)
    )
end

"""
    reverse(fsa)

Return the reversal of fsa.
"""
Base.reverse(fsa::AbstractFSA) = FSA(ω(fsa), T(fsa)', α(fsa), ρ(fsa), λ(fsa))

"""
    renorm(fsa::FSA)

For a state q, multiply each arc's weight by the inverse of the sum of
all the arcs' weight leaving q. In the resulting FSA, sum of the all
the arcs'weight of a state sum up to the semiring one.
"""

function renorm(fsa::FSA{K}) where K
    Z = one(K) ./ (sum(T(fsa), dims=2) .+ ω(fsa))
    FSM(
        inv(sum(α(fsa)) + ρ) * α(fsa),
        T(fsa) .* Z,
        ω(fsa) .* Z[:,1],
        ρ,
        λ(fsa)
    )
end

