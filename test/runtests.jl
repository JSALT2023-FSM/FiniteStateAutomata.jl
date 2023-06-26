using FiniteStateAutomata
import LogExpFunctions: logaddexp
import LinearAlgebra: dot
using Test

@testset "Boolean semiring" begin
    x, y = one(BoolSemiring), zero(BoolSemiring)
    @test ! val(zero(x))
    @test val(one(x))
    @test val(x)
    @test ! val(y)
    @test val(x ⊕ y)
    @test val(x ⊕ x)
    @test ! val(y ⊕ y)
    @test ! val(x ⊗ y)
    @test val(x ⊗ x)
    @test ! val(y ⊗ y)
end

@testset "Logarithmic semiring" begin
    for T in [Float64, Float32]
        a = 2.4
        K = LogSemiring{T,a}
        x, y = K(2), K(3)
        @test val(x ⊕ y) ≈ logaddexp(a*val(x), a*val(y)) / a
        @test val(x ⊗ y) ≈ val(x) + val(y)
        @test val(x ⊘ y) ≈ val(x) - val(y)
        @test val(inv(x)) ≈ -val(x)
        @test val(zero(x)) == T(-Inf)
        @test val(one(x)) == T(0)
    end
end

@testset "Probability semiring" begin
    for T in [Float32, Float64]
        K = ProbSemiring{T}
        x, y = K(2), K(3)
        @test val(x ⊕ y) ≈ val(x) + val(y)
        @test val(x ⊗ y) ≈ val(x) * val(y)
        @test val(x ⊘ y) ≈ val(x) / val(y)
        @test val(inv(x)) ≈ inv(val(x))
        @test val(zero(x)) == T(0)
        @test val(one(x)) == T(1)
    end
end

@testset "Tropical semiring" begin
    for T in [Float32, Float64]
        K = TropicalSemiring{T}
        x, y = K(2), K(3)
        @test val(x ⊕ y) ≈ max(val(x), val(y))
        @test val(x ⊗ y) ≈ val(x) + val(y)
        @test val(x ⊘ y) ≈ val(x) - val(y)
        @test val(inv(x)) ≈ -val(x)
        @test val(zero(x)) == T(-Inf)
        @test val(one(x)) == T(0)
    end
end

