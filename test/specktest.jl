using Random
using CUDA
using FiniteStateAutomata
using SparseArrays

@testset "spECK Sparse matrix -  sparse matrix" begin
    sizes = [512, 1000]
    repeats = 10
    for i in sizes
        for r in 1:repeats
            A = sprand(Float32, i, i, 0.01);
            B = sprand(Float32, i, i, 0.01);

            As = FiniteStateAutomata.sparse(A)
            Bs = FiniteStateAutomata.sparse(B)

            Cx = A * B

            cuAs = to_gpu(As)
            cuBs = to_gpu(Bs)

            cuCq  = Float32.(to_cpu(speck_spgemm(cuAs, cuBs)))
            Cq  = Float32.(speck_spgemm(As, Bs))

            @test Cx ≈ cuCq
            @test Cx ≈ Cq
        end
    end
end
