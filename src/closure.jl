# SPDX-License-Identifier: CECILL-2.1

struct ClosureFST{K,L,TA<:AbstractFST{K,L}} <: AbstractFST{K,L}
    A::TA
    plus::Bool
end

α(cA::ClosureFST) = α(cA.A)
T(cA::ClosureFST) = SPUV(T(cA.A), hcat(ω(cA)), vcat(transpose(α(cA.A))))
ω(cA::ClosureFST) = ω(cA.A)
ρ(cA::ClosureFST{K}) where K = iszero(ρ(cA.A)) && ! cA.plus ? one(K) : ρ(cA.A)
λ(cA::ClosureFST) = λ(cA.A)

closure(A::TransducerOrAcceptor; plus = false) = ClosureFST(A, plus)

