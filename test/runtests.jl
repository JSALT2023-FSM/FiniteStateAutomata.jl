using FiniteStateAutomata
import LogExpFunctions: logaddexp
import LinearAlgebra: dot
using Test

@testset "Semirings" begin
    include("semirings.jl")
end

<<<<<<< HEAD
try
    @testset "Semirings" begin
        include("semirings.jl")
    end
catch e
    @warn "Exception occured in Semirings tests"
=======
@testset "Linear Algebra" begin
    include("matrixtest.jl")
end

@testset "FST" begin
    include("fst.jl")
end

@testset "Reading / Writing FSM format" begin
    include("io.jl")
>>>>>>> draw
end
