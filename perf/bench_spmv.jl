using BenchmarkTools
using FiniteStateAutomata
using SuiteSparseGraphBLAS


bench_mv(dim) = begin
    A = Array(sprand(Float32, dim, dim, 0.05));
    x = Array(sprand(Float32, 1, dim, 0.1));
   
    t = @benchmark $x * $A
    t
end

bench_spmv_default(dim) = begin
    s = sprand(Float32, dim, dim, 0.05);
    v = sprand(Float32, 1, dim, 0.1);

    t = @benchmark $v * $s
    t
end

bench_spmv_gbmatrix(dim) = begin
    s = sprand(Float32, dim, dim, 0.05);
    v = sprand(Float32, 1, dim, 0.1);

    s = GBMatrix(s); v = GBMatrix(v);

    t = @benchmark $v * $s
    t
end

bench_spmv_our(dim) = begin
    s = sprand(Float32, dim, dim, 0.05);
    v = sprand(Float32, 1, dim, 0.1);

    s = sparse(s); v = sparsevec(v);

    t = @benchmark $v * $s
    t
end

io = IOContext(stdout)
sizes_dense  = [100, 200, 500, 1000, 5000, 10000]
sizes_sparse = [100, 200, 500, 1000, 5000, 10000, 100000]

for size in sizes_dense
    println(size, " Julia dense")
    t = bench_mv(size)
    show(io, MIME("text/plain"), t); println("")
end

for size in sizes_sparse
    println(size, " Julia Sparse")
    t = bench_spmv_default(size)
    show(io, MIME("text/plain"), t); println("")
end


for size in sizes_sparse
    println(size, " GBMatrix")
    t = bench_spmv_gbmatrix(size)
    show(io, MIME("text/plain"), t); println("")
end


for size in sizes_sparse
    println(size, " Our own")
    t = bench_spmv_our(size)
    show(io, MIME("text/plain"), t); println("")
end


