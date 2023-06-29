function kronecker_coo(X::SparseMatrixCOO,Y::SparseMatrixCOO)
    rowval = []
    colval = []
    nzval = []
    for (ii, jj, v) in zip(X.rowval, X.colval, X.nzval)
        append!(rowval,ii+Y.rowval-1)
        append!(colval,jj+Y.colval-1)
        append!(nzcal,v*Y)
    end
end