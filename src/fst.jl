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
    T::TransitionMatrix{K}
    ω::AbstractVector{K}
    ρ::K
    λ::AbstractVector{L}
end

FST(A::AbstractFST) = FST(α(A), T(A), ω(A), ρ(A), λ(A))

# Returns a `Q` vector from a list of tuple `(i, v)`.`K` is element
# type of the matrix.
function _spv_from_list(K, Q, tuples)
    I, V = Int[], K[]
    for (i, v) in tuples
        push!(I, i)
        push!(V, v)
    end
    sparsevec(I, V, Q)
end

function FST(K, initweights, transmat, finalweights, statelabels, ϵweight = zero(K))
    Q = length(statelabels)
    FST(
        _spv_from_list(K, Q, initweights),
        transmat,
        _spv_from_list(K, Q, finalweights),
        ϵweight,
        statelabels
    )
end

α(A::AbstractFST) = A.α
T(A::AbstractFST) = A.T
ω(A::AbstractFST) = A.ω
ρ(A::AbstractFST) = A.ρ
λ(A::AbstractFST) = A.λ

