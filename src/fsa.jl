
"""
    struct FSA{K,L} <: AbstractFSA{K,L}
        α::AbstractSparseVector{K}
        T::AbstractSparseMatrix{K}
        ω::AbstractSparseVector{K}
        λ::AbstractVector{L}
    end

Generic Finite State Automaton.
"""
struct FSA{K,L} <: AbstractFSA{K,L}
    α::AbstractSparseVector{K}
    T::AbstractSparseMatrix{K}
    ω::AbstractSparseVector{K}
    λ::AbstractVector{L}
end

