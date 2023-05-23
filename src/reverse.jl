# SPDX-License-Identifier: CECILL-2.1

struct ReversedFST{K,L,TA<:AbstractFST{K,L}} <: AbstractFST{K,L}
    A::TA
end

ReverseFST(rA::ReversedFST) = rA.A

α(rA::ReversedFST) = ω(rA.A)
T(rA::ReversedFST) = copy(T(rA.A) |> transpose)
ω(rA::ReversedFST) = α(rA.A)
ρ(rA::ReversedFST) = ρ(rA.A)
λ(rA::ReversedFST) = λ(rA.A)

"""
    reverse(A)

Reverse all the arcs in A.
"""
Base.reverse(A::AbstractFST) = ReversedFST(A)
