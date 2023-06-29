# SPDX-License-Identifier: CECILL-2.1
# abstract type AbstractSparseMatrixCOO{S} <: AbstractMatrix{S} end

struct SparseMatrixCOO{S} <: AbstractSparseMatrixCOO{S}
    m::Int
    n::Int
    rowval::Vector{Int}
    colval::Vector{Int}
    nzval::Vector{S}
end

# Create an uninitialized sparse matrix with internal buffers to
# store `nzc` non-zero values.
SparseMatrixCOO(S::Type{<:Semiring}, m, n, nzc) = SparseMatrixCOO(
    m,
    n,
    Vector{Int}(undef, nzc),
    Vector{Int}(undef, nzc),
    Vector{S}(undef, nzc)
)

#= Array interface =#

Base.size(X::SparseMatrixCOO) = (X.m, X.n)

function Base.getindex(X::SparseMatrixCOO{S}, i::Int, j::Int) where S
    for (ii, jj, v) in zip(X.rowval, X.colval, X.nzval)
        if ii == i && jj == j
            return v
        end
    end
    return zero(S)    
end

function Base.getindex(X::SparseMatrixCOO{S}, i::Int, ::Colon) where S
    ii = findall(x->x==i, X.rowval)
    transpose(SparseVector(size(X, 2), X.colval[ii], X.nzval[ii]))    
end

#= Sparse array interface =#

rowvals(X::SparseMatrixCOO) = X.rowval
colvals(X::SparseMatrixCOO) = X.colval
nonzeros(X::SparseMatrixCOO) = X.nzval

function sparse_coo(I::AbstractVector, J::AbstractVector, V::AbstractVector{S},
                m = maximum(I), n = maximum(J)) where S
    p = sortperm(I)
    colval = J[p]
    rowval = I[p]
    nzval = V[p]
    SparseMatrixCOO(m, n, rowval, colval, nzval)
end

sparse_coo(I::AbstractVector, J::AbstractVector, val, m = maximum(I), n = maximum(J)) =
    sparse(I, J, repeat([val], length(J)), m, n)

sparse_coo(X::SparseMatrixCSR) = sparse_coo(findnz(X)..., X.m, X.n)

nnz(X::SparseMatrixCOO) = length(X.nzval)

# TODO: complete if needed
# nzrange(X::SparseMatrixCOO, r) = X.rowptr[r]:(X.rowptr[r+1]-1)
# nnz(X::SparseMatrixCOO, r) = X.rowptr[r+1]-X.rowptr[r]

function findnz(X::SparseMatrixCOO)
    X.rowval, X.colval, X.nzval
end

function blockdiag(X::SparseMatrixCOO{S}...) where S
    mX = Int[ size(x, 1) for x in X ]
    nX = Int[ size(x, 2) for x in X ]
    m = sum(mX)
    n = sum(nX)

    rowval = vcat([x.rowval.+(i-1)*x.m for (i,x) in enumerate(X)]...)
    colval = vcat([x.colval.+(i-1)*x.n for (i,x) in enumerate(X) ]...)
    nzval = vcat([x.nzval for (i,x) in enumerate(X) ]...)
    # @show rowval, colval, nzval

    SparseMatrixCOO(m, n, rowval, colval, nzval)
end

