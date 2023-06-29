Base.:*(a::Semiring, X::SparseMatrixCOO) =
    SparseMatrixCOO(size(X)..., collect(rowvals(X)), collect(colvals(X)), a .⊗ collect(nonzeros(X)))
Base.:*(X::SparseMatrixCOO, a::Semiring) =
    SparseMatrixCOO(size(X)..., collect(rowvals(X)), collect(colvals(X)), collect(nonzeros(X)) .⊗ a)
