# SPDX-License-Identifier: CECILL-2.1

using ChainRulesCore
using LinearAlgebra
#using Semirings
using CUDA
using Adapt

export
    # Types
    SparseMatrixCSR,
    SparseMatrixCOO,
    SparseTensorCOO,
    # SparseMatrices,
    sparsematrices,

    # Sparse API
    getrowptr,
    colvals,
    rowvals,
    nonzeros,
    sparsevec,
    sparse_csr,
    sparse_coo,
    sparsetensor_coo,
    nzrange,
    nnz,
    findnz,    
    blockdiag


include("abstracttypes.jl")
include("sparseaccumulator.jl")
include("sparsevector.jl")
include("sparsematrix_csr.jl")
include("sparsematrix_coo.jl")
include("sparsetensor_coo.jl")
include("sparsematrices.jl")
include("kron.jl")
include("linalg.jl")
include("linalg_coo.jl")

#=====================================================================#
# Power series of a matrix.
#=====================================================================#
export powerseries

include("power.jl")


#TODO: Move elsewhere
findparam(ex::Type{<:AbstractArray{T}}) where {T} = T

findbasetype(ex::Type{<:AbstractArray}) = Base.typename(ex).wrapper


