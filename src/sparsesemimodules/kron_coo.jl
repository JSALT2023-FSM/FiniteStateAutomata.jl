
function kron_coo(A,B)
    S = eltype(A)
    ma, mb, na, nb = A.m, B.m, A.n, B.n
    Q = nnz(A)*nnz(B)
    rows = zeros(Int,Q)
    cols = zeros(Int,Q)
    vals = Vector{S}(undef, Q)
    c = 1
    for (i,j,a) in zip(A.rows, A.cols, A.vals)
        for (k,l,b) in zip(B.rows, B.cols, B.vals)
            rows[c] = (i-1)*mb+k
            cols[c] = (j-1)*nb+l
            vals[c] = a*b
            c+=1
        end
    end
    SparseMatrixCOO{S, Int}(ma*mb, na*nb, rows, cols, vals)
end

function kron_coo_mt(A,B)
    S = eltype(A)
    ma, mb, na, nb = A.m, B.m, A.n, B.n
    Q = nnz(A)*nnz(B)
    rows = zeros(Int,Q)
    cols = zeros(Int,Q)
    vals = Vector{S}(undef, Q)
    c = 1

    Threads.@threads for (acoo,bcoo) in collect(Iterators.product(zip(A.rows, A.cols, A.vals), zip(B.rows, B.cols, B.vals)))
        i,j,a = acoo
        k,l,b = bcoo
        rows[c] = (i-1)*mb+k
        cols[c] = (j-1)*nb+l
        vals[c] = a*b
        c+=1
    end
    SparseMatrixCOO{S, Int}(ma*mb, na*nb, rows, cols, vals)
end

function kron_coo_noalloc(A,B)
    rows = Vector{Int}()
    cols = Vector{Int}()
    S = eltype(A)
    vals = Vector{S}()
    ma, mb, na, nb = A.m, B.m, A.n, B.n
    for (i,j,a) in zip(A.rows, A.cols, A.vals)
        for (k,l,b) in zip(B.rows, B.cols, B.vals)
            push!(rows,(i-1)*mb+k)
            push!(cols,(j-1)*nb+l)
            push!(vals, a*b)
        end
    end
    SparseMatrixCOO{S, Int}(ma*mb, na*nb, rows, cols, vals)
end

function kron_coo2dict(A,B)  
    D = Dict()
    ma, mb, na, nb = A.m, B.m, A.n, B.n 
    for (i,j,a) in zip(A.rows, A.cols, A.vals)
        for (k,l,b) in zip(B.rows, B.cols, B.vals)
            c = a*b
            D[((i-1)*mb+k,(j-1)*nb+l )] = c
        end
    end
    D
end

Base.kron(A::AbstractSparseMatrixCOO{S}, B::AbstractSparseMatrixCOO{S}) where S =
    kron_coo(A, B)
