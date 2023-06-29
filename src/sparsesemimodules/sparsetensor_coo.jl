# SPDX-License-Identifier: CECILL-2.1

struct SparseTensorCOO{S} <: AbstractSparseTensorCOO{S}
    size::Tuple{Int}
    coo::Matrix{Int}
    nzval::Vector{S}
end

# Create an uninitialized sparse matrix with internal buffers to
# store `nzc` non-zero values.
SparseTensorCOO(S::Type{<:Semiring}, size, nzc) = SparseTensorCOO(
    size,
    Matrix{Int}(undef, nzc),
    Vector{S}(undef, nzc)
)

#= Array interface =#

Base.size(X::SparseTensorCOO) = X.size

# function Base.getindex(X::SparseTensorCOO{S}, i::Int, j::Int) where S
#     for (ii, jj, v) in zip(X.rowval, X.colval, X.nzval)
#         if ii == i && jj == j
#             return v
#         end
#     end
#     return zero(S)    
# end

# function Base.getindex(X::SparseTensorCOO{S}, i::Int, ::Colon) where S
#     ii = findall(x->x==i, X.rowval)
#     transpose(SparseVector(size(X, 2), X.colval[ii], X.nzval[ii]))    
# end

#= Sparse array interface =#

# rowvals(X::SparseTensorCOO) = X.rowval
# colvals(X::SparseTensorCOO) = X.colval
# nonzeros(X::SparseTensorCOO) = X.nzval

# function sparsetensor_coo(coo::AbstractMatrix, V::AbstractVector{S}, shape) where S
#     SparseTensorCOO(rowval, colval, nzval)
# end

nnz(X::SparseTensorCOO) = length(X.nzval)

# TODO: complete if needed
# nzrange(X::SparseTensorCOO, r) = X.rowptr[r]:(X.rowptr[r+1]-1)
# nnz(X::SparseTensorCOO, r) = X.rowptr[r+1]-X.rowptr[r]

# function findnz(X::SparseTensorCOO)
#     X.rowval, X.colval, X.nzval
# end

# function blockdiag(X::SparseTensorCOO{S}...) where S
#     mX = Int[ size(x, 1) for x in X ]
#     nX = Int[ size(x, 2) for x in X ]
#     m = sum(mX)
#     n = sum(nX)

#     rowval = vcat([x.rowval.+(i-1)*x.m for (i,x) in enumerate(X)]...)
#     colval = vcat([x.colval.+(i-1)*x.n for (i,x) in enumerate(X) ]...)
#     nzval = vcat([x.nzval for (i,x) in enumerate(X) ]...)
#     # @show rowval, colval, nzval

#     SparseTensorCOO(m, n, rowval, colval, nzval)
# end

