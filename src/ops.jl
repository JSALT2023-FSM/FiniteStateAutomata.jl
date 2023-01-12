# SPDX-License-Identifier: MIT

"""
    union(fsa1[, fsa2, ...])
    fsa1 ∪ fsa2

Return the union of the given FSA.
"""
function Base.union(fsa1::AbstractFSA{K,L}, fsa2::AbstractFSA{K,L}) where {K,L}
    FSA(
        vcat(fsa1.α, fsa2.α),
        blockdiag(fsa1.T, fsa2.T),
        vcat(fsa1.ω, fsa2.ω),
        fsa1.ρ + fsa2.ρ,
        vcat(fsa2.λ, fsa2.λ)
    )
end
Base.union(fsa1::AbstractFSA{K,L}, fsaN::AbstractFSA{K,L}...) where {K,L} =
    foldl(union, fsaN, init = fsa1)

"""
    cat(fsa1[, fsa2, ...])

Return the concatenation of the given FSA.
"""
function Base.cat(fsa1::AbstractFSA{K,L}, fsa2::AbstractFSA{K,L}) where {K,L}
    D1, D2 = size(fsa1.T, 1), size(fsa2.T, 1)
    FSA(
        vcat(fsa1.α, zero(K) * fsa2.α),
        [fsa1.T (fsa1.ω * fsa2.α');
         spzeros(K, D2, D1) fsa2.T],
        vcat(zero(K) * fsa1.ω, fsa2.ω),
        fsa1.ρ * fsa2.ρ,
        vcat(fsa2.λ, fsa2.λ)
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
        fsa.α,
        fsa.T + fsa.ω * fsa.α',
        fsa.ω,
        iszero(fsa.ρ) && ! plus ? one(K) : fsa.ρ,
        fsa.λ
    )
end

"""
    reverse(fsa)

Return the reversal of fsa.
"""
Base.reverse(fsa::AbstractFSA) = FSA(fsa.ω, fsa.T', fsa.α, fsa.ρ, fsa.λ)
