# SPDX-License-Identifier: CECILL-2.1

struct CatFST{K,L,TA<:AbstractFST{K,L},TB<:AbstractFST{K,L}} <: AbstractFST{K,L}
    A::TA
    B::TB
end

α(cA::CatFST) = vcat(α(cA.A), ρ(cA.A) * α(cA.B))
T(cA::CatFST{K}) where K = [
    T(cA.A) (ω(cA.A) * transpose(α(cA.B)));
    spzeros(K, nstates(cA.B), nstates(cA.A)) T(cA.B)
]

ω(cA::CatFST) = vcat(ρ(cA.B) * ω(cA.A), ω(cA.B))
ρ(cA::CatFST) = ρ(cA.A) ⊗ ρ(cA.B)
λ(cA::CatFST) = vcat(λ(cA.A), λ(cA.B))

"""
    cat(A1[, A2, ...])

Return the concatenation of the given FSTs.
"""
Base.cat(A::TransducerOrAcceptor, B::TransducerOrAcceptor) = CatFST(A, B)
Base.cat(A1::AbstractFST{K}, AN::AbstractFST{K}...) where K = foldl(cat, AN, init = A1)

