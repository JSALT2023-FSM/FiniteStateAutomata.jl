# SPDX-License-Identifier: CECILL-2.1

"""
    relabel(fn, fst)

Create a new FST where each labelis map each label
"""

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
