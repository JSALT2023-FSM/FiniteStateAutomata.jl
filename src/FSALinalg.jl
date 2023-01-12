# SPDX

module FSALinalg

using LinearAlgebra
using SparseArrays

export
    # FSA types
    FSA,

    # FSA operations
    # Base.union
    # Base.cat
    closure

include("abstractfsa.jl")
include("fsa.jl")
include("ops.jl")

end
