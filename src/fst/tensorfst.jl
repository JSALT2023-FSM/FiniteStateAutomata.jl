# SPDX-License-Identifier: CECILL-2.1

"""
    struct TensorFST{S} <: AbstractExpandedFST{S}
        M::AbstractVector{<:AbstractMatrix{S}}
        α::AbstractVector{S}
        ω::AbstractVector{S}
    end

Finite State Transducer using sparse arrays for storage.
"""
struct TensorFST{S} <: AbstractExpandedFST{S}
    M::AbstractArray{S, 4}
    α::AbstractVector{S}
    ω::AbstractVector{S}
end

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


M(fst::TensorFST) = fst.M
α(fst::TensorFST) = fst.α
ω(fst::TensorFST) = fst.ω

numstates(fst::TensorFST) = size(fst.M, 1)

numarcs(fst::TensorFST) = sum(!iszero, fst.M)

states(fst::TensorFST) = 1:numstates(fst)

arcs(fst::TensorFST, q) = [(coo[1], coo[2], coo[3], fst.M[q,coo]) for coo in findall(!iszero, fst.M[q, :, :, :])]

initstate(fst::TensorFST) = findall(!iszero, fst.α)[1]

finalweight(fst::AbstractFST, q) =  fst.ω[q]

