# SPDX-License-Identifier: CECILL-2.1

# using ChainRulesCore
using LinearAlgebra
#using Semirings
using CUDA
using Adapt

export
    # Types
    SparseMatrixCOO,
    hasitem,
    # SparseMatrixCSR,
    # SparseMatrices,

    # Sparse API
    # getrowptr,
    # colvals,
    # rowvals,
    # nonzeros,
    # sparsevec,
    # sparse,
    # nzrange,
    # nnz,
    # findnz,
    # blockdiag,
    tocoo,
    todense


include("abstracttypes.jl")
include("sparseaccumulator.jl")
include("sparsevector.jl")
include("sparsematrix_coo.jl")
# include("sparsematrix_csr.jl")
# include("sparsematrices.jl")
include("kron_coo.jl")
include("conversions.jl")
# include("linalg.jl")

#=====================================================================#
# Power series of a matrix.
#=====================================================================#
export powerseries

include("power.jl")


#TODO: Move elsewhere
findparam(ex::Type{<:AbstractArray{T}}) where {T} = T

findbasetype(ex::Type{<:AbstractArray}) = Base.typename(ex).wrapper


