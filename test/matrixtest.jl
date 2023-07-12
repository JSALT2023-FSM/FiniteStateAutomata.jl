u:qa
sing Random
using SparseArrays
using CUDA

@testset "CUDA Sparse matrix and vector" begin
    #sizes = [1, 3, 5, 50, 512, 5000, 50000]
    sizes = [512]
    repeats = 1
    @show "CUDA Sparse matrix and vector"
    for i in sizes
        for r in 1:repeats
            A = sprand(Float32, i, i, 0.01);
            x = sprand(Float32, 1, i, 0.1);

            xs = FiniteStateAutomata.sparsevec(x)
            As = FiniteStateAutomata.sparse(A)
            ys = xs * As

            cuxs = to_gpu(xs)
            cuas = to_gpu(As)
            cuys = cuxs * cuas
            CUDA.@time  cuxs * cuas


            o = to_cpu(cuys)
            @test vec(o) ≈ vec(ys)
      end
  end
end

@testset "CUDA Sparse matrix and matrix" begin
    rng = Random.MersenneTwister(1234)
    #sizes = [1, 3, 5, 50, 500, 5000]
    sizes = [5000]
    repeats = 1
    @show "CUDA Sparse matrix and matrix"
    for i in sizes
        for r in 1:repeats
            A = sprand(Float32, i, i, 0.01);
            B = sprand(Float32, i, i, 0.01);
            
            @show i, r
            CUDA.@time w = A*B
            As = FiniteStateAutomata.sparse(A)
            Bs = FiniteStateAutomata.sparse(B)
            #@time Ys = FiniteStateAutomata.sparse(As * Bs)

            cubs = to_gpu(Bs)
            cuas = to_gpu(As)
            cuys = cuas * cubs
            CUDA.@time cuas * cubs

            o= to_cpu(cuys)
            CUDA.@time p = A * B
            #@show o
            #@show p
            @test vec(o) ≈ vec(p)
            #@test vec(o) ≈ vec(Ys)
      end
  end
end

@testset "Sparse vector-matrix multiplication" begin
    sizes = [1, 3, 5, 50, 500, 5000]
    repeats = 10
    for i in sizes
        for r in 1:repeats
            A = Float32.(rand([-1.0, 2,0, 0, 0, 0],(i,i)))
            x = Float32.(rand([-1.0, 0], (1, i)))
            y = x * A

            xs = FiniteStateAutomata.sparsevec(x)
            As = FiniteStateAutomata.sparse(A)
            ys = xs * As

            @test vec(y) ≈ vec(y)
        end
    end
end

@testset "Mostly empty matrix" begin
  A = Float32[ 1 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
  As = FiniteStateAutomata.sparse(A)

  y2 = SparseMatrixCSR(As.m, As.n, As.rowptr, As.colval, As.nzval)
  @test y2 == A
end

@testset "Sparse to float32 and back" begin
  A = Float32[ 1 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
  As = FiniteStateAutomata.sparse(A)

  A2 = Float32.(As)

  @test A2 ≈ A
end

