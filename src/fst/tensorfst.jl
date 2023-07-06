# SPDX-License-Identifier: CECILL-2.1

"""
    struct TensorFST{S,T,P<:Tuple{Int,Int,Int,Int}} <: ExpandedFST{S}
        M::T
        α::AbstractVector{S}
        ω::AbstractVector{S}
    end

Finite State Transducer using sparse arrays for storage.
`T` is type of 4D tensor, e.g. `Array{S, 4}` for dense FST.
`M` is a 4D tensor encoding the arcs. The first dimension encodes
the source state, the second dimension encodes the destination state,
the third dimension encodes input label and the fourth dimension is the
output label. The the type is parameterize by a tuple indicating the
actual "orientation" of the tensor. `α` encodes the arcs leaving the
super initial state `0` and `ω` encodes the final weights.
"""
struct TensorFST{S,T<:AbstractArray{S,4},P} <: ExpandedFST{S}
    M::T
    α::AbstractVector{S}
    ω::AbstractVector{S}

    function TensorFST(M::AbstractArray{S,N}, α, ω, p) where {S,N}
        N == 4 || throw(ArgumentError("Tensor has to be 4-dimensional."))
        length(α) == size(M, 1) || throw(ArgumentError("α has different #states than M"))
        length(α) == length(ω) || throw(ArgumentError("ω has different #states than α"))
        sum(findall(!iszero, α)) == 1 || throw(ArgumentError("TensorFST has to have only 1 initstate"))
        new{S,typeof(M),p}(M, α, ω)
    end
end

# If no permutation is given assumes default permutation.
TensorFST(M, α, ω) = TensorFST(M, α, ω, (1, 2, 3, 4))

#=====================================================================#
# Implementation of the ExpandedFST interface.
#=====================================================================#

numstates(fst::TensorFST) = size(fst.M, 1)
numarcs(fst::TensorFST) = sum(!iszero, fst.M)
states(fst::TensorFST) = 1:numstates(fst)
arcs(fst::TensorFST, q) = [
    (coo[1], coo[2], coo[3], fst.M[q,coo])
    for coo in findall(!iszero, fst.M[q, :, :, :])
]
initstate(fst::TensorFST) = findall(!iszero, fst.α)[1]
finalweight(fst::TensorFST, q) = fst.ω[q]

#=====================================================================#
# TensorFST interface.
#=====================================================================#

"""
    M(fst)

 A four-dimensional tensor representing the arcs (source node, destination node, and weight)
"""
M(fst::TensorFST) = fst.M


"""
    α(fst)

Vector of the initial weights.
"""
α(fst::TensorFST) = fst.α


"""
    ω(fst)

Vector of the final weights.
"""
ω(fst::TensorFST) = fst.ω

