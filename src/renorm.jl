# SPDX-License-Identifier: CECILL-2.1

struct RenormFST{K,L,TA<:AbstractFST{K,L}} <: AbstractFST{K,L}
    A::TA
end

function α(rA::RenormFST)
    Z = inv(sum(nonzeros(α(rA.A))) ⊕ ρ(rA.A))
    Z * α(rA.A)
end

function T(rA::RenormFST)
    Z = inv.(sum(eachcol(T(rA.A))) + ω(rA.A))
    Z .⊗ T(rA.A)
end

function ω(rA::RenormFST)
    Z = inv.(sum(eachcol(T(rA.A))) + ω(rA.A))
    Z[:, 1] .* ω(rA.A)
end

function ρ(rA::RenormFST)
    Z = inv(sum(nonzeros(α(rA.A))) ⊕ ρ(rA.A))
    Z * ρ(rA.A)
end

λ(rA::RenormFST) = λ(rA.A)

"""
    renorm(A)

Local renormalization. For a state q, multiply each arc's weight by the
inverse of the sum of all the arcs' weight leaving q. In the resulting
FST, sum of the all the arcs'weight of a state sum up to the semiring
one.
"""
renorm(A::TransducerOrAcceptor) = RenormFST(A)

