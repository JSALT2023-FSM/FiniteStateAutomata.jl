# SPDX-License-Identifier: CECILL-2.1

"""
    abstract type AbstractSparseVector{S} <: AbstractVector{S}

Abstract base type for sparse vectors.
"""
abstract type AbstractSparseVector{S} <: AbstractVector{S} end

"""
    abstract type AbstractSparseMatrixCSR{S} <: AbstractMatrix{S}

Abstract base type for sparse matrices in [CSR](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)) format.
"""
abstract type AbstractSparseMatrixCSR{S} <: AbstractMatrix{S} end


"""
    getrowptr(M)

Return the row pointer of the matrix.
"""
getrowptr(::AbstractSparseMatrixCSR)

"""
    colvals(M::AbstractSparseMatrixCSR)

Return the column indices of the non-zero values of `M`.
"""
colvals

"""
    rowvals(x::AbstractSparseVector)

Return the row indices of the non-zero values of `x`.
"""
rowvals

"""
    nonzeros(M)

Return the non-zero values.
"""
nonzeros(::AbstractSparseMatrixCSR)

"""
    sparsevec(I, V, nrows)

Create a sparse vector.
"""
sparsevec

"""
    sparse_csr(I, J, V, nrows, ncols)

Create a sparse CSR matrix.
"""
sparse_csr



"""
    nzrange(x)

Return an iterator over the non-zero values of `x`.

    nzrange(M, r)

Return an iterator over the non-zero values of the row `r`.
"""
nzrange(::AbstractSparseVector)

"""
    nzrange(M, r)

Return an iterator over the non-zero values of the row `r`.
"""
nzrange(::AbstractSparseMatrixCSR, r)

"""
    I, V = findnz(x)
    I, J, V = findnz(M)

Convert the sparse vector/matrix in a coordinate format.
"""
findnz

"""
    nnz(x)

Return the number of non-zero values in `x`.
"""
nnz

"""
    abstract type AbstractSparseMatrixCOO{S} <: AbstractMatrix{S}

Abstract base type for sparse matrices in [COO](https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)) format.
"""
abstract type AbstractSparseMatrixCOO{S} <: AbstractMatrix{S} end

"""
    sparse_coo(I, J, V, nrows, ncols)

Create a sparse COO matrix.
"""
sparse_coo



"""
    abstract type AbstractSparseMatrixCOO{S} <: AbstractMatrix{S}

Abstract base type for sparse matrices in [COO](https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)) format.
"""
abstract type AbstractSparseTensorCOO{S,N} <: AbstractArray{S,N} end
