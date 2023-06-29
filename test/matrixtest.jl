### A Pluto.jl notebook ###
# v0.19.26



@testset "Sparse vector-matrix multiplication" begin
    sizes = [1, 3, 5, 50, 5000]
    repeats = 10
    for i in sizes
        for r in 1:repeats
            A = Float32.(rand([-1.0, 2,0, 0, 0, 0],(i,i)))
            x = Float32.(rand([-1.0, 0], (1, i)))

            y = x * A
        
        
            xs = sparsevec(x)
            As = sparse(A)

            #println(xs)
            #println(As)
            #println(x)
            #println(A)

            ys = xs * As
            #println(typeof(As))
            #@show y, typeof(y)
            #@show ys, typeof(ys)
            @test ys ≈ y
        end
    end
end

@testset "Mostly empty matrix" begin
  A = Float32[ 1 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
  As = sparse(A)

  y2 = SparseMatrixCSR(As.m, As.n, As.rowptr, As.colval, As.nzval)
  @test y2 == A
end

@testset "Sparse to float32 and back" begin
  A = Float32[ 1 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
  As = sparse(A)
  
  A2 = Float32.(As)

  @test A2 ≈ A
end

