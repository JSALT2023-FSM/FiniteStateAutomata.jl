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

f(L, A, w, i) = begin
    l = i <= nstates(A) ? λ(A)[i] : one(L)
    S = UnionConcatSemiring{L}
    if iszero(w)
        return ProductSemiring((w, zero(S)))
    else
        return ProductSemiring((w, S(Set(L[l]))))
    end
end

for K in Ks, L in Ls
    @testset verbose=true "FSA - $K" begin
        α1 = sparsevec([1, 2], K[2, 3], 3)
        T1 = sparse([1, 1, 2], [2, 3, 3], K[1, 2, 3], 3, 3)
        ω1 = sparsevec([3], K(5), 3)
        ρ1 = K(0.6)
        λ1 = [L("a"), L("b"), L("c")]
        A1 = FSA(α1, T1, ω1, ρ1, λ1)
        B1 = convert((w, l) -> f(L, A1, w, l), A1)

        α2 = sparsevec([1, 3], K[3, 5], 4)
        T2 = sparse([1, 1, 2, 3], [2, 3, 4, 4], K[4, 3, 2, 2], 4, 4)
        ω2 = sparsevec([2, 4], K[5, 6], 4)
        ρ2 = zero(K)
        λ2 = [L("a"), L("b"), L("c"), L("d")]
        A2 = AcyclicFSA(α2, T2, ω2, ρ2, λ2)
        B2 = convert((w, l) -> f(L, A2, w, l), A2)

        A3 = FSA(
            sparsevec([1, 4], one(K), 6),
            sparse([1, 2, 4, 5], [2, 3, 5, 6], one(K), 6, 6),
            sparsevec([3, 6], one(K), 6),
            one(K),
            [L("a"), L("b"), L("c"), L("a"), L("e"), L("c")]
        )
        B3 = convert((w, l) -> f(L, A3, w, l), A3)

        A4 = FSA(
            sparsevec([1], one(K) * 2, 4),
            sparse([1, 1, 2, 3], [2, 3, 4, 4], one(K), 4, 4),
            sparsevec([4], one(K) * 2, 4),
            one(K),
            [L("a"), L("b"), L("e"), L("c")]
        )
        B4 = convert((w, l) -> f(L, A4, w, l), A4)

        Aϵ = FSA(spzeros(K, 0), spzeros(K, 0, 0), spzeros(K, 0), one(K), L[])
        Bϵ = convert((w, l) -> f(L, Aϵ, w, l), Aϵ)

        @testset verbose=true "properties" begin
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
        end

        # Number of iteration to sum over the FSA.
        n = max(nstates(A1), nstates(A2))

        @testset verbose=true "union" begin
            A12 = union(A1, A2)
            B12 = convert((w, l) -> f(L, A12, w, l), A12)
            cs_B12 = sum(B12; n)
            cs_B1_B2 = sum(B1; n) + sum(B2; n)
            @test cs_B12.tval[1] ≈ cs_B1_B2.tval[1]
            @test cs_B12.tval[2] == cs_B1_B2.tval[2]
            @test typeof(union(A2, A2)) <: AbstractAcyclicFSA
        end

        @testset verbose=true "concatenation" begin
            A12 = cat(A1, A2)
            B12 = convert((w, l) -> f(L, A12, w, l), A12)
            cs_B12 = sum(B12; n = 2n)
            cs_B1_B2 = sum(B2; n) * sum(B1; n)
            @test cs_B12.tval[2] ⊇ cs_B1_B2.tval[2]
            @test typeof(cat(A2, A2)) <: AbstractAcyclicFSA
        end

        # The number of iteration for the cumulative sum depends on the
        # structure of the FSA for the test to pass.
        # It should be 3 times the maximum path length of B2.
        n = 6

        @testset verbose=true "closure" begin
            B2p_n3 = B2 ∪ cat(B2, B2) ∪ cat(B2, B2, B2)
            B2_n3 = B2 ∪ cat(B2, B2) ∪ cat(B2, B2, B2) ∪ Bϵ
            cB2p = closure(B2; plus = true)
            cB2 = closure(B2; plus = false)
            cs_B2p_n3 = sum(B2p_n3; n)
            cs_B2_n3 = sum(B2_n3; n)
            cs_cB2p = sum(cB2p; n)
            cs_cB2 = sum(cB2; n)
            @test cs_cB2p.tval[1] ≈ cs_B2p_n3.tval[1]
            @test cs_cB2p.tval[2] == cs_B2p_n3.tval[2]
            @test cs_cB2.tval[1] ≈ cs_B2_n3.tval[1]
            @test cs_cB2.tval[2] == cs_B2_n3.tval[2]
        end

        @testset verbose=true "reversal" begin
            rA2 = A2 |> reverse
            rB2 = convert((w, l) -> f(L, rA2, w, l), rA2)
            s1 = val(sum(B2; n).tval[2])
            s2 = Set((StringMonoid ∘ reverse ∘ val).((val ∘ sum)(rB2)[2]))
            @test s1 == s2
            @test typeof(rA2) <: AbstractAcyclicFSA
        end

        @testset verbose=true "renorm" begin
            @test Base.isapprox(val(sum(A1 |> renorm; n = 100)), val(one(K)), atol=1e-6)
            @test Base.isapprox(val(sum(A2 |> renorm; n = 100)), val(one(K)), atol=1e-6)
        end

        @testset verbose=true "(co)accessible" begin
            A = FSA(
                sparsevec([1, 2], one(K), 5),
                sparse([2, 3], [3, 3], one(K), 5, 5),
                sparsevec([3, 4], one(K), 5),
                zero(K),
                [L("$i") for i in 1:5]
            )
            @test all(accessible(A) .== [true, true, true, false, false])
            @test all(coaccessible(A) .== [false, true, true, true, false])
        end

        @testset verbose=true "globalrenorm" begin
            A = AcyclicFSA(A2) |> globalrenorm
            @test Base.isapprox(val(sum(A; n = 100)), val(one(K)), atol=1e-6)
            @test typeof(A) <: AbstractAcyclicFSA
        end

        @testset verbose=true "minimize" begin
            mA3 = A3 |> minimize
            mB3 = convert((w, l) -> f(L, mA3, w, l), mA3)
            s1 = sum(B4; n = 4)
            s2 = sum(mB3; n = 4)
            @test s1.tval[1] ≈ s2.tval[1]
            @test s1.tval[2] == s2.tval[2]
            @test typeof(AcyclicFSA(A3) |> minimize) <: AbstractAcyclicFSA
        end
    end
end

