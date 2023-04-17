# SPDX-License-Identifier: CECILL-2.1

struct StateMappedFST{K,L,TA<:AbstractFST{K,L}} <: AbstractFST{K,L}
    A::TA
    M::AbstractMatrix{K}
    λ::AbstractVector{L}
end

statemap(A::AbstractFST, M, λ) = StateMappedFST(A, M, λ)
statemap(mA::StateMappedFST, M, λ) = StateMappedFST(mA.A, M' * M, λ)

function statemap(f, A::AbstractFST{K}, λ) where K
    I, J = Int[], Int[]
    for q in 1:nstates(A)
        J_q = f(q)
        I = vcat(I, repeat([q], length(J_q)))
        J = vcat(J, J_q)
    end
    println(I)
    M = sparse(I, J, one(K), nstates(A), maximum(J))
    statemap(A, M, λ)
end

function Base.filter(f, A::AbstractFST{K}) where K
    I = findall(f, 1:nstates(A))
    M = sparse(I, 1:length(I), one(K), nstates(A), length(I))
    statemap(A, M, λ(A)[I])
end

function ChainRulesCore.rrule(::typeof(statemap), A::AbstractFST{K}, M::AbstractMatrix{K}, l) where K
    Y = StateMappedFST(A, M, l)
    pullback(ΔY) = (NoTangent(), StateMappedFST(ΔY, M', λ(A)), NoTangent(), NoTangent())
    Y, pullback
end

α(fA::StateMappedFST) = fA.M' * α(fA.A)
T(fA::StateMappedFST) = fA.M' * T(fA.A) * fA.M
ω(fA::StateMappedFST) = fA.M' * ω(fA.A)
ρ(fA::StateMappedFST) = ρ(fA.A)

