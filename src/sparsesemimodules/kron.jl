# SPDX-License-Identifier: CECILL-2.1

struct KronMatrixCSR{S} <: AbstractSparseMatrixCSR{S}
    A::AbstractSparseMatrixCSR{S}
    B::AbstractSparseMatrixCSR{S}
end

Base.kron(A::AbstractSparseMatrixCSR{S}, B::AbstractSparseMatrixCSR{S}) where S =
    KronMatrixCSR(A, B)

#= Array interface =#

Base.size(M::KronMatrixCSR) = (size(M.A, 1) * size(M.B, 1), size(M.A, 2) * size(M.B, 2))

function Base.getindex(M::KronMatrixCSR, i::Integer, j::Integer)
    rA = (i-1) ÷ size(M.B, 1) + 1
    rB = (i-1) % size(M.B, 1) + 1
    cA = (j-1) ÷ size(M.B, 2) + 1
    cB = (j-1) % size(M.B, 2) + 1
    M.A[rA,cA] ⊗ M.B[rB,cB]
end

#= Sparse array API =#

struct KronRowptrIterator
    M::KronMatrixCSR
end

Base.length(it::KronRowptrIterator) = size(it.M, 1) + 1
Base.iterate(it::KronRowptrIterator, state=1) =
    state > length(it) ? nothing : (it[state], state+1)

function Base.getindex(it::KronRowptrIterator, i::Integer)
    i == size(it.M, 1) + 1 && return nnz(it.M) + 1
    rA = (i-1) ÷ size(it.M.B, 1) + 1
    rB = (i-1) % size(it.M.B, 1) + 1
    ptrA = getrowptr(it.M.A)
    ptrB = getrowptr(it.M.B)
    nzb_rA = ptrA[rA]-1
    nzb_rB = ptrB[rB]-1
    nzc_rA = ptrA[rA+1]-ptrA[rA]
    nzb_rA * nnz(it.M.B) + nzc_rA * nzb_rB + 1
end

getrowptr(M::KronMatrixCSR) = KronRowptrIterator(M)

struct KronColIterator
    M::KronMatrixCSR
end

function Base.getindex(it::KronColIterator, t)
    cA = colvals(it.M.A)[first(t)]
    cB = colvals(it.M.B)[last(t)]
    (cA - 1) * size(it.M.A, 2) + cB
end

Base.length(it::KronColIterator) = nnz(it.M.A) * nnz(it.M.B)

Base.iterate(it::KronColIterator, state=1) =
    state > length(it) ? nothing : (it[state], state+1)

colvals(M::KronMatrixCSR) = KronColIterator(M)

struct KronNzvalIterator
    nzA
    nzB
end

function Base.getindex(it::KronNzvalIterator, t::Tuple)
    vA = it.nzA[first(t)]
    vB = it.nzB[last(t)]
    vA ⊗ vB
end

Base.length(it::KronNzvalIterator) = length(it.vA) * length(it.nzB)
Base.iterate(it::KronNzvalIterator, state=1) =
    state > length(it) ? nothing : (it[state], state+1)

nonzeros(M::KronMatrixCSR) = KronNzvalIterator(nonzeros(M.A), nonzeros(M.B))

struct KronNzrangeIterator
    rangeA
    rangeB
end

Base.length(it::KronNzrangeIterator) = length(it.rangeA) * length(it.rangeB)

function Base.iterate(it::KronNzrangeIterator, state=(first(it.rangeA),first(it.rangeB)))
    i, j = state
    (i > last(it.rangeA)) && return nothing
    j > last(it.rangeB)  && return iterate(it, (i+1, first(it.rangeB)))
    (state, (i, j+1))
end

function nzrange(M::KronMatrixCSR, r::Integer)
    rA = (r-1) ÷ size(M.B, 1) + 1
    rB = (r-1) % size(M.B, 1) + 1
    ptrA = getrowptr(M.A)
    ptrB = getrowptr(M.B)
    KronNzrangeIterator(ptrA[rA]:(ptrA[rA+1]-1), ptrB[rB]:(ptrB[rB+1]-1))
end

nnz(M::KronMatrixCSR) = nnz(M.A) * nnz(M.B)

