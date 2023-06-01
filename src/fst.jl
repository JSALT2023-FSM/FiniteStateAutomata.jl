# SPDX-License-Identifier: CECILL-2.1

# Returns a `Q` x `P` sparse matrix from a list of arcs.
# An arc is defined by a triplet `(i, j, v)`. `K` is element type
# of the matrix.
function _spm_from_list(K, Q, P, arcs)
    I, J, V = Int[], Int[], K[]
    for (src, dest, weight) in arcs
        push!(I, src)
        push!(J, dest)
        push!(V, weight)
    end
    sparse(I, J, V, Q, P)
end

function TransitionMatrix(K, Q, direct_arcs, in_factors, factors, out_factors)
    S = _spm_from_list(K, Q, Q, direct_arcs)

    # Number of factors.
    P = max(
        maximum(t -> t[2], in_factors),
        maximum(t -> max(t[1], t[2]), factors),
        maximum(t -> t[1], out_factors)
    )

    U = _spm_from_list(K, Q, P, in_factors)
    E = _spm_from_list(K, P, P, factors)
    V = _spm_from_list(K, P, Q, out_factors)

    TransitionMatrix(S, U, E, V)
end

"""
    struct FST{K,L} <: AbstractWFST{K,L}
        α::AbstractSparseVector{K}
        T::AbstractSparseMatrix{K}
        ω::AbstractSparseVector{K}
        ρ::K
        λ::AbstractVector{L}
    end

Generic Finite State Automaton.
"""
struct WFST{K,L} <: AbstractWFST{K,L}
    α::AbstractVector{K}
    T::TransitionMatrix{K}
    ω::AbstractVector{K}
    ρ::K
    λ::AbstractVector{L}
end

FST(A::AbstractWFST) = FST(α(A), T(A), ω(A), ρ(A), λ(A))

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

function WFST(; semiring, initweights, arcs, finalweights, statelabels,
             ϵweight = zero(semiring), infactors = [], factors = [], outfactors = [])
    Q = length(statelabels)

    # Number of factors.
    P = max(
        maximum(t -> t[2], infactors; init = 0),
        maximum(t -> max(t[1], t[2]), factors; init = 0),
        maximum(t -> t[1], outfactors; init = 0)
    )

    S = _spm_from_list(semiring, Q, Q, arcs)
    U = _spm_from_list(semiring, Q, P, infactors)
    E = _spm_from_list(semiring, P, P, factors)
    V = _spm_from_list(semiring, P, Q, outfactors)

    WFST(
        _spv_from_list(semiring, Q, initweights),
        TransitionMatrix(S, U, E,V),
        _spv_from_list(semiring, Q, finalweights),
        ϵweight,
        statelabels
    )
end

α(A::AbstractWFST) = A.α
T(A::AbstractWFST) = A.T
ω(A::AbstractWFST) = A.ω
ρ(A::AbstractWFST) = A.ρ
λ(A::AbstractWFST) = A.λ

