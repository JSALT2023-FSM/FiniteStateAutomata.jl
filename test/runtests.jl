using FiniteStateAutomata
import LogExpFunctions: logaddexp
import LinearAlgebra: dot
using Test

@testset "FST" begin
    include("fst.jl")
end

@testset "Semirings" begin
    include("semirings.jl")
end
