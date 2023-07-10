
function kron_coo(A,B)
    rows = Vector{Int}()
    cols = Vector{Int}()
    S = eltype(A)
    vals = Vector{S}()
    ma, mb, na, nb = A.m, B.m, A.n, B.n 
    for (i,j,a) in zip(A.rows, A.cols, A.vals)
        for (k,l,b) in zip(B.rows, B.cols, B.vals)
            c = a*b
            push!(rows,(i-1)*mb+k)
            push!(cols,(j-1)*nb+l)
            push!(vals, c)
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
