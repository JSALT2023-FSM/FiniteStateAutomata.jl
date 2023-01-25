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

α(fsa::AbstractFSA) = parent(fsa).α
T(fsa::AbstractFSA) = parent(fsa).T
ω(fsa::AbstractFSA) = parent(fsa).ω
ρ(fsa::AbstractFSA) = parent(fsa).ρ
λ(fsa::AbstractFSA) = parent(fsa).λ

function Base.convert(f::Function, A::AbstractFSA{K,L}) where {K,L}
    l = λ(A)
    ρ = f(emptystring(A), one(L))
    U = typeof(ρ)
    α = sparsevec(
        initstates(A)[1],
        [f(v, l[i]) for (i, v) in zip(initstates(A)...)],
        nstates(A)
    )
    T = sparse(
        edges(A)[1],
        edges(A)[2],
        [f(v, l[j]) for (i, j, v) in zip(edges(A)...)],
        nstates(A),
        nstates(A)
    )
    ω = sparsevec(
        finalstates(A)[1],
        [f(v, one(L)) for (i, v) in zip(finalstates(A)...)],
        nstates(A)
    )

    FSA{U,L}(α, T, ω, ρ, l)
end

struct AcyclicFSA{K,L} <: AbstractAcyclicFSA{K,L}
    fsa::FSA{K,L}
end

AcyclicFSA(α, T, ω, ρ, λ) = AcyclicFSA(FSA(α, T, ω, ρ, λ))

Base.parent(A::AcyclicFSA) = A.fsa

Base.convert(f, A::AbstractAcyclicFSA) = AcyclicFSA(convert(f, A.fsa))

