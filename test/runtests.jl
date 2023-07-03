using FiniteStateAutomata
import LogExpFunctions: logaddexp
import LinearAlgebra: dot
using Test

@testset "Semirings" begin
    include("semirings.jl")
end

@testset "Linear Algebra" begin
    include("matrixtest.jl")
end

@testset "FST" begin
    include("fst.jl")
end

