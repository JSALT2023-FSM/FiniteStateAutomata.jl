# SPDX-License-Identifier: CECILL-2.1

"""
    struct TensorFST{S, T} <: ExpandedFST{S}
        M::T
        α::AbstractVector{S}
        ω::AbstractVector{S}
    end

Finite State Transducer using sparse arrays for storage.
`T` is type of 4D tensor, e.g. `Array{S, 4}` for dense FST.
"""
struct TensorFST{S, T<:AbstractArray{S, 4}} <: ExpandedFST{S}
    M::T
    α::AbstractVector{S}
    ω::AbstractVector{S}

    function TensorFST(M::AbstractArray{S, N}, α, ω) where {S, N}
        N == 4 || throw(ArgumentError("Matrix has to be 4-dimensional."))
        length(α) == size(M, 1) || throw(ArgumentError("α has differet #states than M"))
        length(α) == length(ω) || throw(ArgumentError("ω has differet #states than α"))
        sum(findall(!iszero, α)) == 1 || throw(ArgumentError("TensorFST has to have only 1 initstate"))
        new{S, typeof(M)}(M, α, ω)
    end
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

arcs(fst::TensorFST, q) = [
    (coo[1], coo[2], coo[3], fst.M[q,coo])
    for coo in findall(!iszero, fst.M[q, :, :, :])
]

initstate(fst::TensorFST) = findall(!iszero, fst.α)[1]

finalweight(fst::AbstractFST, q) =  fst.ω[q]

