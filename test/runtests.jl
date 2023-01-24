using FSALinalg
using Semirings
using SparseArrays
using Test

Ks = [
    LogSemiring{Float64},
    LogSemiring{Float32},
    ProbSemiring{Float64},
    LogSemiring{Float32},
    TropicalSemiring{Float64},
    TropicalSemiring{Float32}
]

Ls = [
    StringMonoid,
]

f(L, w, l) = begin
    S = UnionConcatSemiring{L}
    if iszero(w)
        return ProductSemiring((w, zero(S)))
    else
        return ProductSemiring((w, S(Set(L[l]))))
    end
end

@testset "FSA" begin
    for K in Ks, L in Ls
        α1 = sparsevec([1, 2], K[2, 3], 3)
        T1 = sparse([1, 1, 2], [2, 3, 3], K[1, 2, 3], 3, 3)
        ω1 = sparsevec([3], K(5), 3)
        ρ1 = K(0.6)
        λ1 = [L("a"), L("b"), L("c")]
        A1 = FSA(α1, T1, ω1, ρ1, λ1)
        B1 = convert((w, l) -> f(L, w, l), A1)

        α2 = sparsevec([1], K[3], 4)
        T2 = sparse([1, 1, 2, 3], [2, 3, 4, 4], K[4, 3, 2, 1], 4, 4)
        ω2 = sparsevec([2, 4], K[5, 6], 4)
        ρ2 = zero(K)
        λ2 = [L("a"), L("b"), L("c"), L("d")]
        A2 = FSA(α2, T2, ω2, ρ2, λ2)
        B2 = convert((w, l) -> f(L, w, l), A2)

        Aϵ = FSA(spzeros(K, 0), spzeros(K, 0, 0), spzeros(K, 0), one(K), L[])
        Bϵ = convert((w, l) -> f(L, w, l), Aϵ)

        @test all(α(A1) .≈ α1)
        @test all(T(A1) .≈ T1)
        @test all(ω(A1) .≈ ω1)
        @test all(ρ(A1) ≈ ρ1)
        @test all(λ(A1) .== λ1)

        @test all(α(A2) .≈ α2)
        @test all(T(A2) .≈ T2)
        @test all(ω(A2) .≈ ω2)
        @test all(ρ(A2) ≈ ρ2)
        @test all(λ(A2) .== λ2)

        # union
        B12 = convert((w, l) -> f(L, w, l), union(A1, A2))
        cs_B12 = cumsum(B12)
        cs_B1_B2 = cumsum(B1) + cumsum(B2)
        @test cs_B12.tval[1] ≈ cs_B1_B2.tval[1]
        @test cs_B12.tval[2] == cs_B1_B2.tval[2]

        # concatenation
        B12 = convert((w, l) -> f(L, w, l), cat(A1, A2))
        cs_B12 = cumsum(B12)
        cs_B1_B2 = cumsum(B2) * cumsum(B1)
        @test cs_B12.tval[1] ≈ cs_B1_B2.tval[1]
        @test cs_B12.tval[2] == cs_B1_B2.tval[2]

        # The number of iteration for the cumulative sum depends on the
        # structure of the FSA for the test to pass.
        # It should be 3 times the maximum path length of B2.
        n = 6

        # closure
        B2p_n3 = B2 ∪ cat(B2, B2) ∪ cat(B2, B2, B2)
        B2_n3 = B2 ∪ cat(B2, B2) ∪ cat(B2, B2, B2) ∪ Bϵ
        cB2p = closure(B2; plus = true)
        cB2 = closure(B2; plus = false)
        cs_B2p_n3 = cumsum(B2p_n3; n = n)
        cs_B2_n3 = cumsum(B2_n3; n = n)
        cs_cB2p = cumsum(cB2p; n = n)
        cs_cB2 = cumsum(cB2; n = n)
        @test cs_cB2p.tval[1] ≈ cs_B2p_n3.tval[1]
        @test cs_cB2p.tval[2] == cs_B2p_n3.tval[2]
        @test cs_cB2.tval[1] ≈ cs_B2_n3.tval[1]
        @test cs_cB2.tval[2] == cs_B2_n3.tval[2]
    end
end

