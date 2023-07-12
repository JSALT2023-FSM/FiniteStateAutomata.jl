using FiniteStateAutomata

lib="libjulia_specklib.so"

i = 10
A = Float32.(rand([-1.0, 2,0, 0, 0, 0],(i,i)))
x = Float32.(rand([-1.0, 0], (1, i)))
y = x * A

xs = FiniteStateAutomata.sparsevec(x)
As = FiniteStateAutomata.sparse(A)
Cx = A * A
Cs = FiniteStateAutomata.sparse(Cx)

function speck_spgemm(A::SparseMatrixCSR{Float32}, B::SparseMatrixCSR{Float32})
    nnzA = nnz(A)
    nnzB = nnz(B)

    width = Ref{Csize_t}(0)
    C_rowptr = Ref{Ptr{Cuint}}(0)
    C_colval = Ref{Ptr{Cuint}}(0)
    C_nzval = Ref{Ptr{Float32}}(0)

    Arowptr = UInt32.(A.rowptr .-1)
    Acolval = UInt32.(A.colval .-1)

    Browptr = UInt32.(B.rowptr .-1)
    Bcolval = UInt32.(B.colval .-1)

    @ccall lib.julia_multiply_float_cpu(A.m::Csize_t, A.n::Csize_t, B.n::Csize_t,
                                   nnzA::Csize_t,
                                   Arowptr::Ptr{Cuint},
                                   Acolval::Ptr{Cuint},
                                   A.nzval::Ptr{Cfloat},
                                   nnzB::Csize_t,
                                   Browptr::Ptr{Cuint},
                                   Bcolval::Ptr{Cuint},
                                   B.nzval::Ptr{Cfloat},
                                   width::Ref{Csize_t},
                                   C_rowptr::Ptr{Ptr{Cuint}},
                                   C_colval::Ptr{Ptr{Cuint}},
                                   C_nzval::Ptr{Ptr{Float32}})::Cvoid

    nnz0 = Int(width[])
    rowptr = Int.(unsafe_wrap(Array, C_rowptr[], A.m+1; own=false)) .+ 1
    colval = Int.(unsafe_wrap(Array, C_colval[], nnz0; own=false)) .+ 1
    nzval = copy(unsafe_wrap(Array, C_nzval[], nnz0; own=false))

    @ccall lib.julia_float_free_cpumem(C_rowptr[]::Ptr{Cuint},
                                       C_colval[]::Ptr{Cuint},
                                       C_nzval[]::Ptr{Float32})::Cvoid
    return SparseMatrixCSR(A.m, B.n, rowptr, colval, nzval)
end

@show nnz(As)
@show As.rowptr
@show As.colval
@show As.nzval

@show nnz(Cs)
@show Cs.rowptr
@show Cs.colval
@show Cs.nzval

Cq  = Float32.(speck_spgemm(As, As))

@show Cq

@show Cx
@show Cx .≈ Cq
@show Cx ≈ Cq
