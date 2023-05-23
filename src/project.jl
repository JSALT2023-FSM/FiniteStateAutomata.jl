# SPDX-License-Identifier: CECILL-2.1

struct ProjectedFST{K,L,TA<:Transducer{K}} <: AbstractFST{K,L}
    A::TA
    λ::AbstractVector{L}
end

α(pA::ProjectedFST) = α(pA.A)
T(pA::ProjectedFST) = T(pA.A)
ω(pA::ProjectedFST) = ω(pA.A)
ρ(pA::ProjectedFST) = ρ(pA.A)
λ(pA::ProjectedFST) = pA.λ


_project(op, A::Acceptor) = A
_project(op, A::Transducer) = ProjectedFST(A, [op(l) for l in λ(A)])
Π₁(A::TransducerOrAcceptor) = _project(first, A)
Π₂(A::TransducerOrAcceptor) = _project(last, A)

