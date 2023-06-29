# SPDX-License-Identifier: CECILL-2.1

"""
    struct SparseFST{S,L} <: AbstractFST{S,L}
        M::AbstractVector{<:AbstractMatrix{S}}
        α::AbstractVector{S}
        ω::AbstractVector{S}
        λ::AbstractVector{S}
    end

Finite State Transducer using sparse arrays for storage.
"""
struct SparseFST{S,L} <: AbstractFST{S,L}
    M::SparseMatrices{<:SparseMatrixCOO{S}}
    α::SparseVector{S}
    ω::SparseVector{S}
    λ::AbstractVector{L}
end

M(fst::SparseFST) = fst.M
α(fst::SparseFST) = fst.α
ω(fst::SparseFST) = fst.ω
λ(fst::SparseFST) = fst.λ

function arcs(fst::SparseFST)
    retval = []
    for (l, Mᵢ) in enumerate(fst.M)        
        for (i, j, v) in zip(findnz(Mᵢ)...)
           push!(retval, (i, j, fst.λ[l], v))
        end
    end
    retval
end

narcs(fst::SparseFST) = nnz(parent(fst.M))

states(fst::SparseFST) = [(i, fst.α[i], fst.ω[i]) for i in 1:nstates(fst)]


