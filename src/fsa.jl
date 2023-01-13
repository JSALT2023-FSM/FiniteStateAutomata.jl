# SPDX-License-Identifier: CECILL-2.1

"""
    struct FSA{K,L} <: AbstractFSA{K,L}
        α::AbstractSparseVector{K}
        T::AbstractSparseMatrix{K}
        ω::AbstractSparseVector{K}
        ρ::K
        λ::AbstractVector{L}
    end

Generic Finite State Automaton.
"""
struct FSA{K,L} <: AbstractFSA{K,L}
    α::AbstractVector{K}
    T::AbstractMatrix{K}
    ω::AbstractVector{K}
    ρ::K
    λ::AbstractVector{L}
end

α(fsa::FSA) = fsa.α
T(fsa::FSA) = fsa.T
ω(fsa::FSA) = fsa.ω
ρ(fsa::FSA) = fsa.ρ
λ(fsa::FSA) = fsa.λ
