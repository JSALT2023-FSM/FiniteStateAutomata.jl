# SPDX-License-Identifier: CECILL-2.1

struct ReverseFST{K,L,TA<:AbstractFST{K,L}} <: AbstractFST{K,L}
    A::TA
end

ReverseFST(rA::ReverseFST) = rA.A

α(rA::ReverseFST) = ω(rA.A)
T(rA::ReverseFST) = T(rA.A)'
ω(rA::ReverseFST) = α(rA.A)
ρ(rA::ReverseFST) = ρ(rA.A)
λ(rA::ReverseFST) = λ(rA.A)

"""
    reverse(A)

Reverse all the arcs in A.
"""
Base.reverse(A::AbstractFST) = ReverseFST(A)
