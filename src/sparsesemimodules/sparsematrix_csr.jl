# SPDX-License-Identifier: CECILL-2.1

struct SparseMatrixCSR{S} <: AbstractSparseMatrixCSR{S}
    m::Int
    n::Int
    rowptr::Vector{Int}
    colval::Vector{Int}
    nzval::Vector{S}
end

# Create an uninitialized sparse matrix with internal buffers to
# store `nzc` non-zero values.
SparseMatrixCSR(S::Type{<:Semiring}, m, n, nzc) = SparseMatrixCSR(
    m,
    n,
    Vector{Int}(undef, m+1),
    Vector{Int}(undef, nzc),
    Vector{S}(undef, nzc)
)

function gather!(X::SparseMatrixCSR, acc::SparseAccumulator, r::Integer, nzcur::Integer)
    # [TODO: check dimensions] && throw(DimensionMismatch())

    X.rowptr[r] = nzcur
    X.rowptr[r+1] = nzcur + nnz(acc)
    for n in 1:nnz(acc)
        i = acc.ind[n]
        X.colval[nzcur + n - 1] = i
        X.nzval[nzcur + n - 1] = acc.val[i]
    end

    X
end

#= Array interface =#

Base.size(X::SparseMatrixCSR) = (X.m, X.n)

function Base.getindex(X::SparseMatrixCSR{S}, i::Int, j::Int) where S
    ii = nzrange(X, i)
    isempty(ii) && return zero(S)
    nnzy = colvals(X)[ii]
    fj = findfirst(a -> a==j, nnzy)
    (isnothing(fj)) ? zero(S) : nonzeros(X)[ii[1]-1+fj]
end

function Base.getindex(X::SparseMatrixCSR{S}, i::Int, ::Colon) where S
    ii = nzrange(X, i)
    transpose(SparseVector(size(X, 2), colvals(X)[ii], nonzeros(X)[ii]))
end

#= Sparse array interface =#

getrowptr(X::SparseMatrixCSR) = X.rowptr
colvals(X::SparseMatrixCSR) = X.colval
nonzeros(X::SparseMatrixCSR) = X.nzval

function sparse(I::AbstractVector, J::AbstractVector, V::AbstractVector{S},
                m = maximum(I), n = maximum(J)) where S
    p = sortperm(I)
    colval = J[p]

    counts = fill!(similar(I, m+1), 0)
    counts[1] = 1
    for i in I
        counts[i+1] += 1
    end
    rowptr = cumsum(counts)

    SparseMatrixCSR(m, n, rowptr, colval, V[p])
end

sparse(I::AbstractVector, J::AbstractVector, val, m = maximum(I), n = maximum(J)) =
    sparse(I, J, repeat([val], length(J)), m, n)

# TODO: commenting this out for now, as it is not used anywhere
# sparse(X::SparseMatrixCSR) = sparsecsr(findnz(X)..., X.m, X.n)

nzrange(X::SparseMatrixCSR, r) = X.rowptr[r]:(X.rowptr[r+1]-1)
nnz(X::SparseMatrixCSR) = length(X.nzval)
nnz(X::SparseMatrixCSR, r) = X.rowptr[r+1]-X.rowptr[r]

function findnz(X::SparseMatrixCSR)
    I = similar(X.colval, 0)
    for i in 1:size(X.rowptr, 1)-1
        I = vcat(I, repeat([i], nnz(X, i)))
    end
    I, X.colval, X.nzval
end

function blockdiag(X::SparseMatrixCSR{S}...) where S
    num = length(X)
    mX = Int[ size(x, 1) for x in X ]
    nX = Int[ size(x, 2) for x in X ]
    m = sum(mX)
    n = sum(nX)

    rowptr = Vector{Int}(undef, m+1)
    nnzX = Int[ nnz(x) for x in X ]
    nnz_res = sum(nnzX)
    colval = Vector{Int}(undef, nnz_res)
    nzval = Vector{S}(undef, nnz_res)

    nnz_sofar = 0
    nX_sofar = 0
    mX_sofar = 0
    for i = 1 : num
        rowptr[(1 : mX[i] + 1) .+ mX_sofar] = getrowptr(X[i]) .+ nnz_sofar
        colval[(1 : nnzX[i]) .+ nnz_sofar] = colvals(X[i]) .+ nX_sofar
        nzval[(1 : nnzX[i]) .+ nnz_sofar] = nonzeros(X[i])
        nnz_sofar += nnzX[i]
        nX_sofar += nX[i]
        mX_sofar += mX[i]
    end
    rowptr[m+1] = nnz_sofar + 1

    SparseMatrixCSR(m, n, rowptr, colval, nzval)
end

