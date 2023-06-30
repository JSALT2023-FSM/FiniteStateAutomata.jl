# SPDX-License-Identifier: CECILL-2.1

#= left / right scalar multiplication =#

Base.:*(a::Semiring, x::SparseVector) = SparseVector(x.n, x.nzind, a .⊗ x.nzval)
Base.:*(x::SparseVector, a::Semiring) = SparseVector(x.n, x.nzind, x.nzval .⊗ a)
Base.:*(a::Semiring, X::SparseMatrixCSR) =
    SparseMatrixCSR(size(X)..., collect(getrowptr(X)), collect(colvals(X)), a .⊗ collect(nonzeros(X)))
Base.:*(X::SparseMatrixCSR, a::Semiring) =
    SparseMatrixCSR(size(X)..., collect(getrowptr(X)), collect(colvals(X)), collect(nonzeros(X)) .⊗ a)

#= Element-wise addition =#

function Base.:+(x::SparseVector{S}, y::SparseVector{S}) where S
    size(x, 1) != size(y, 1) && throw(DimensionMismatch())

    N = length(x)
    acc = SparseAccumulator(S, N, nnz(x) + nnz(y))

    for i in 1:nnz(x) scatter!(acc, x.nzval[i], x.nzind[i]) end
    for i in 1:nnz(y) scatter!(acc, y.nzval[i], y.nzind[i]) end

    gather!(SparseVector(S, N, nnz(acc)), acc)
end

function Base.:+(X::SparseMatrixCSR{S}, Y::SparseMatrixCSR{S}) where S
    size(X) != size(Y) && throw(DimensionMismatch())

    Xnzval = nonzeros(X)
    Xnzcol = colvals(X)
    Ynzval = nonzeros(Y)
    Ynzcol = colvals(Y)

    Z = SparseMatrixCSR(S , size(X, 1), size(X, 2), nnz(X) + nnz(Y))
    acc = SparseAccumulator(S, size(X, 2), nnz(X) + nnz(Y))
    nzcur = 1
    for r in 1:size(X,1)
        for i in nzrange(X, r) scatter!(acc, Xnzval[i], Xnzcol[i]) end
        for i in nzrange(Y, r) scatter!(acc, Ynzval[i], Ynzcol[i]) end
        gather!(Z, acc, r, nzcur)
        nzcur += nnz(acc)
        reset!(acc)
    end

    nzc = nzcur - 1
    SparseMatrixCSR(Z.m, Z.n, Z.rowptr, Z.colval[1:nzc], Z.nzval[1:nzc])
end

#= dot product =#

function LinearAlgebra.dot(x::SparseVector{S}, y::SparseVector{S}) where S
    size(x, 1) != size(y, 1) && throw(DimensionMismatch())

    acc = SparseAccumulator(S, length(x), nnz(x))
    for i in 1:nnz(x) scatter!(acc, x.nzval[i], x.nzind[i]) end

    I, V = findnz(y)
    s = zero(S)
    for i in 1:nnz(y)
        j = I[i]
        if acc.mask[j]
            s = s ⊕ (acc.val[j] ⊗ V[i])
        end
    end

    s
end

##=== Matrix-vector product ===#

function Base.:*(xᵀ::Transpose{S,<:SparseVector}, A::SparseMatrixCSR{S}) where S
    size(xᵀ, 2) != size(A, 1) && throw(DimensionMismatch())

    x = parent(xᵀ)
    Anzcol = colvals(A)
    Anzval = nonzeros(A)
    xnzcol = rowvals(x)
    xnzval = nonzeros(x)

    acc = SparseAccumulator(S, size(A, 2), size(A, 2))

    for r in nzrange(x)
        for i in nzrange(A, xnzcol[r])
            scatter!(acc, xnzval[r] ⊗ Anzval[i], Anzcol[i])
        end
    end
    y = SparseVector(S, size(A, 2), nnz(acc))
    transpose(gather!(y, acc))

end

