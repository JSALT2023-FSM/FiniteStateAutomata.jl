# SPDX-License-Identifier: CECILL-2.1

"""
    struct FST{K,L} <: AbstractFST{K,L}
        α::AbstractSparseVector{K}
        T::AbstractSparseMatrix{K}
        ω::AbstractSparseVector{K}
        ρ::K
        λ::AbstractVector{L}
    end

Generic Finite State Automaton.
"""
struct FST{K,L} <: AbstractFST{K,L}
    α::AbstractVector{K}
    T::AbstractMatrix{K}
    ω::AbstractVector{K}
    ρ::K
    λ::AbstractVector{L}
end

FST(A::AbstractFST) = FST(α(A), T(A), ω(A), ρ(A), λ(A))

function FST(K, arcs, finalweights, statelabels, ϵweight = zero(K))
    I_α, V_α = Int[], K[]
    I_T, J_T, V_T = Int[], Int[], K[]
    for (src, dest, weight) in arcs
        if src == 0
            push!(I_α, dest)
            push!(V_α, weight)
        else
            push!(I_T, src)
            push!(J_T, dest)
            push!(V_T, weight)
        end
    end

    I_ω, V_ω = Int[], K[]
    for (state, weight) in finalweights
        push!(I_ω, state)
        push!(V_ω, weight)
    end

    Q = length(statelabels)
    FST(
        sparsevec(I_α, V_α, Q),
        sparse(I_T, J_T, V_T, Q, Q),
        sparsevec(I_ω, V_ω, Q),
        ϵweight,
        statelabels
    )
end

α(A::AbstractFST) = A.α
T(A::AbstractFST) = A.T
ω(A::AbstractFST) = A.ω
ρ(A::AbstractFST) = A.ρ
λ(A::AbstractFST) = A.λ

