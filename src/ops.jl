# SPDX-License-Identifier: CECILL-2.1

"""
    W(A)

Return the total weight of `A`, i.e. the ``\\oplus``-sum of all the
path's weight in `A`.
"""
W(A::AbstractWFST) =
    ρ(A) ⊕ (transpose(α(A)) * MatrixPowerSum(T(A)) * ω(A)) # TODO optimize with dot(., ., .)

function closure(A::AbstractWFST; plus = false)
    K = eltype(α(A))
    TA = T(A)
    TB = TransitionMatrix(
        TA.S,
        hcat(TA.U, ω(A)),
        blockdiag(TA.E, spzeros(K, 1, 1)),
        vcat(TA.V, transpose(α(A)))
    )
    ρB = iszero(ρ(A)) && ! plus ? one(K) : ρ(A)
    WFST(α(A), TB, ω(A), ρB, λ(A))
end

"""
    cat(A1[, A2, ...])

Return the concatenation of the given FSTs.
"""
function Base.cat(A::AbstractWFST, B::AbstractWFST)
    K = promote_type(eltype(α(A)), eltype(α(B)))
    TA, TB = T(A), T(B)
    TC = TransitionMatrix(
        blockdiag(TA.S, TB.S),
        hcat(blockdiag(TA.U, TB.U), vcat(ω(A), spzeros(K, nstates(B)))),
        blockdiag(blockdiag(TA.E, TB.E), spzeros(K, 1, 1)),
        vcat(blockdiag(TA.V, TB.V), transpose(vcat(spzeros(K, nstates(B)), α(B)))),
    )

    WFST(
        vcat(α(A), ρ(A) * α(B)),
        TC,
        vcat(ω(A) * ρ(B), ω(B)),
        ρ(A) ⊗ ρ(B),
        vcat(λ(A), λ(B))
    )
end
Base.cat(A1::AbstractWFST{K}, AN::AbstractWFST{K}...) where K = foldl(cat, AN, init = A1)

Π₁(A::Acceptor) = WFST(α(A), T(A), ω(A), ρ(A), λ(A))
Π₂(A::Acceptor) = WFST(α(A), T(A), ω(A), ρ(A), λ(A))
Π₁(A::AbstractWFST) = WFST(α(A), T(A), ω(A), ρ(A), first.(λ(A)))
Π₂(A::AbstractWFST) = WFST(α(A), T(A), ω(A), ρ(A), last.(λ(A)))

"""
    union(A1[, A2, ...])
    A1 ∪ A2

Return the union of the given FST.
"""
function Base.union(A::AbstractWFST, B::AbstractWFST)
    K = promote_type(eltype(α(A)), eltype(α(B)))
    TA, TB = T(A), T(B)
    TC = TransitionMatrix(
        blockdiag(TA.S, TB.S),
        blockdiag(TA.U, TB.U),
        blockdiag(TA.E, TB.E),
        blockdiag(TA.V, TB.V)
    )

    WFST(vcat(α(A), α(B)), TC, vcat(ω(A), ω(B)), ρ(A) ⊕ ρ(B), vcat(λ(A), λ(B)))
end

Base.union(A1::AbstractWFST{K}, AN::AbstractWFST{K}...) where K = foldl(union, AN, init = A1)
