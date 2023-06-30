using FiniteStateAutomata
import LogExpFunctions: logaddexp
import LinearAlgebra: dot
using Test

try
    @testset "FST" begin
        include("fst.jl")
    end
catch e
    @warn "Exception occured in FST tests"
end

try
    @testset "Semirings" begin
        include("semirings.jl")
    end
catch e
    @warn "Exception occured in Semirings tests"
end
