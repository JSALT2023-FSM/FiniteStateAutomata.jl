# SPDX-License-Identifier: CECILL-2.1

struct SparseMatrices{T<:AbstractMatrix} <: AbstractVector{T}
    rowptr::Vector{Int}
    colptr::Vector{Int}
    matrix::T
    matrixes::Vector{T}
end

sparsematrices(X::AbstractMatrix...) = SparseMatrices(
    cumsum(vcat([1], [size(x, 1) for x in X])),
    cumsum(vcat([1], [size(x, 2) for x in X])),
    blockdiag(X...),
    [x for x in X]
)

Base.parent(X::SparseMatrices) = X.matrix
Base.size(X::SparseMatrices) = (length(X.rowptr) - 1,)

function Base.getindex(X::SparseMatrices{SparseMatrixCSR}, i::Integer)
    rows = X.matrix.rowptr[X.rowptr[i]]:(X.matrix.rowptr[X.rowptr[i+1]]-1)
    m = X.rowptr[i+1] - X.rowptr[i]
    n = X.colptr[i+1] - X.colptr[i]

    SparseMatrixCSR(
        m,
        n,
        X.matrix.rowptr[X.rowptr[i]:X.rowptr[i+1]] .- (X.matrix.rowptr[X.rowptr[i]] - 1),
        X.matrix.colval[rows] .- (X.colptr[i] - 1),
        X.matrix.nzval[rows]
    )
end

function Base.getindex(X::SparseMatrices{S}, i::Integer) where S<:SparseMatrixCOO
    X.matrixes[i]
end

#=====================================================================#
# Broadcasting
#=====================================================================#

Base.BroadcastStyle(::Type{<:SparseMatrices}) = Broadcast.ArrayStyle{SparseMatrices}()

function Base.Broadcast.broadcasted(::typeof(*), a::S, X::SparseMatrices) where S
    SparseMatrices(X.rowptr, X.colptr, a * X.matrix)
end

function Base.Broadcast.broadcasted(::typeof(*), X::SparseMatrices, a::S) where S
    SparseMatrices(X.rowptr, X.colptr, X.matrix * a)
end

function Base.broadcasted(::typeof(*), a::Transpose{S,<:SparseVector},
                                 X::SparseMatrices{S}) where S
    @show "not implemented"
end


