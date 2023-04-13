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

α(A::AbstractFST) = parent(A).α
T(A::AbstractFST) = parent(A).T
ω(A::AbstractFST) = parent(A).ω
ρ(A::AbstractFST) = parent(A).ρ
λ(A::AbstractFST) = parent(A).λ

function Base.convert(f::Function, A::AbstractFST{K,L}) where {K,L}
    ρ = f(emptystring(A), nstates(A) + 1)
    U = typeof(ρ)
    α = sparsevec(
        initstates(A)[1],
        [f(v, i) for (i, v) in zip(initstates(A)...)],
        nstates(A)
    )
    T = sparse(
        edges(A)[1],
        edges(A)[2],
        [f(v, j) for (i, j, v) in zip(edges(A)...)],
        nstates(A),
        nstates(A)
    )
    ω = sparsevec(
        finalstates(A)[1],
        [f(v, nstates(A) + 1) for (i, v) in zip(finalstates(A)...)],
        nstates(A)
    )

    FST{U,L}(α, T, ω, ρ, λ(A))
end

struct AcyclicFST{K,L} <: AbstractAcyclicFST{K,L}
    fsa::AbstractFST{K,L}
end

AcyclicFST(α, T, ω, ρ, λ) = AcyclicFST(FST(α, T, ω, ρ, λ))

Base.parent(A::AcyclicFST) = parent(A.fsa)

Base.convert(f::Function, A::AbstractAcyclicFST{K,L}) where {K,L} =
    AcyclicFST(convert(f, parent(A)))

