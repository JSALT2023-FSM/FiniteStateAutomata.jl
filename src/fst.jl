# SPDX-License-Identifier: CECILL-2.1

"""
    struct FST{K,L} <: AbstractFST{K,L}
        α::AbstractVector{K}
        S::AbstractMatrix{K}
        D::AbstractMatrix{K}
        ω::AbstractSparseVector{K}
        λ::AbstractVector{L}
        isymbols::Dict
        osymbols::Dict
    end

Generic Finite State Automaton.
"""
struct FST{K,L} <: AbstractFST{K,L}
    α::AbstractVector{K}
    S::AbstractMatrix{K}
    D::AbstractMatrix{K}
    ω::AbstractVector{K}
    λ::AbstractVector{L}
end

FST(A::AbstractFST; isymbols = isymbols(A), osymbols = osymbols(A)) =
    FST(α(A), S(A), D(A), ω(A), λ(A), isymbols, osymbols)

α(A::FST) = A.α
S(A::FST) = A.S
D(A::FST) = A.D
ω(A::FST) = A.ω
ρ(A::FST) = A.ρ
λ(A::FST) = A.λ

