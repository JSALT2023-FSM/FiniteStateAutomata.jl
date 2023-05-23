# SPDX-License-Identifier: CECILL-2.1

struct UnionFST{K,L,TA<:AbstractFST{K,L},TB<:AbstractFST{K,L}} <: AbstractFST{K,L}
    A::TA
    B::TB
end

α(uA::UnionFST) = vcat(α(uA.A), α(uA.B))
T(uA::UnionFST) = blockdiag(T(uA.A), T(uA.B))
ω(uA::UnionFST) = vcat(ω(uA.A), ω(uA.B))
ρ(uA::UnionFST) = ρ(uA.A) ⊕ ρ(uA.B)
λ(uA::UnionFST) = vcat(λ(uA.A), λ(uA.B))

"""
    union(A1[, A2, ...])
    A1 ∪ A2

Return the union of the given FST.
"""
Base.union(A::TransducerOrAcceptor, B::TransducerOrAcceptor) = UnionFST(A, B)
Base.union(A1::AbstractFST{K}, AN::AbstractFST{K}...) where K = foldl(union, AN, init = A1)

