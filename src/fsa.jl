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

function Base.convert(f::Function, A::AbstractFSA{K,L}) where {K,L}
    l = λ(A)
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

    FSA(α, T, ω, iszero(ρ(A)) ? zero(eltype(α)) : f(ρ(A), one(L)), l)
end

