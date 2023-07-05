# SPDX-License-Identifier: CECILL-2.1

using LinearAlgebra
#using Semirings
using CUDA
using Adapt

export
    # Types
    SparseMatrixCSR,
    SparseMatrices,

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

#=====================================================================#
# Power series of a matrix.
#=====================================================================#
export powerseries

include("power.jl")


#TODO: Move elsewhere
findparam(ex::Type{<:AbstractArray{T}}) where {T} = T

findbasetype(ex::Type{<:AbstractArray}) = Base.typename(ex).wrapper


