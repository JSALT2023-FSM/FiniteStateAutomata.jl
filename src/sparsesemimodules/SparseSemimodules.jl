# SPDX-License-Identifier: CECILL-2.1

using ChainRulesCore
using LinearAlgebra
#using Semirings
using CUDA
using Adapt

export
    # Types
    SparseMatrixCSR,
    SparseMatrices,

    # CUDA accelerated interface
    CuSparseMatrixCSR,
    CuSparseVectorX,

    # Sparse API
    getrowptr,
    colvals,
    rowvals,
    nonzeros,
    sparsevec,
    sparse,
    nzrange,
    nnz,
    findnz,
    blockdiag



include("abstracttypes.jl")
include("sparseaccumulator.jl")
include("sparsevector.jl")
include("sparsematrix_csr.jl")
include("sparsematrices.jl")
include("kron.jl")
include("linalg.jl")
include("cu_linalg.jl")

export
    to_gpu,
    to_cpu,
    mult_spvspm!   # sparse vector/sparse matrix multiplication
                   # buffer (dense output vector must be provided)
    mult_spmspm!

#=====================================================================#
# Power series of a matrix.
#=====================================================================#
export powerseries

include("power.jl")

include("speck.jl")
export speck_spgemm

#TODO: Move elsewhere
findparam(ex::Type{<:AbstractArray{T}}) where {T} = T

findbasetype(ex::Type{<:AbstractArray}) = Base.typename(ex).wrapper


