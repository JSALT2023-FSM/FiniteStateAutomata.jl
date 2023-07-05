import  SparseArrays.SparseMatrixCSC

function csc2coo(A)
	colptr = A.colptr
	rowval = A.rowval
	vals = A.nzval
	rows = Vector{Int}()
	cols = Vector{Int}()
	for i in 1:length(colptr)-1
		for j in colptr[i]:colptr[i+1]-1
			push!(rows,rowval[j])
			push!(cols,i)
		end
	end
	SparseMatrixCOO{eltype(vals),Int}(A.m,A.n,rows, cols, vals)
end

tocoo(A::SparseMatrixCSC) = csc2coo(A)

function coo2dense(A) 
    D = zeros(eltype(A),(A.m,A.n))
    for k = 1:nnz(A)
        i, j, v = A.rows[k], A.cols[k], A.vals[k]
        D[i,j] = v
    end
    D
end

todense(A::SparseMatrixCOO) = coo2dense(A)