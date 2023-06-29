# SPDX-License-Identifier: CECILL-2.1

"""
    M(fst)

 A four-dimensional tensor representing the arcs (source node, destination node, and weight)     
"""
M(fst::TensorFST)


"""
    α(fst)

Vector of the initial weights.
"""
α(fst::TensorFST)


"""
    ω(fst)

Vector of the final weights.
"""
ω(fst::TensorFST)

"""
    struct TensorFST{S} <: AbstractFST{S}
        M::AbstractVector{<:AbstractMatrix{S}}
        α::AbstractVector{S}
        ω::AbstractVector{S}
    end

Finite State Transducer using sparse arrays for storage.
"""
struct TensorFST{S} <: AbstractFST{S}
    M::SparseMatrices{<:SparseMatrixCSR{S}}
    α::SparseVector{S}
    ω::SparseVector{S}
end

M(fst::TensorFST) = fst.M
α(fst::TensorFST) = fst.α
ω(fst::TensorFST) = fst.ω


function arcs(fst::TensorFST)
    retval = []
    for (l, Mᵢ) in enumerate(fst.M)
        for (i, j, v) in zip(findnz(Mᵢ)...)
           push!(retval, (i, j, fst.λ[l], v))
        end
    end
    retval
end

narcs(fst::TensorFST) = nnz(parent(fst.M))

states(fst::TensorFST) = [(i, fst.α[i], fst.ω[i]) for i in 1:nstates(fst)]

