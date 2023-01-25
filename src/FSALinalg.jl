# SPDX-License-Identifier: CECILL-2.1

module FSALinalg

using LinearAlgebra
using SparseArrays

export
    # Abstract types
    AbstractFSA,
    AbstractAcyclicFSA,

    # concrete types
    AcyclicFSA,
    FSA,

    # Accessors / properties
    α,
    T,
    ω,
    ρ,
    λ,
    nstates,
    nedges,
    initstates,
    edges,
    finalstates,

    # FSA operations
    closure,
    globalrenorm,
    renorm


include("abstractfsa.jl")
include("fsa.jl")
include("ops.jl")

end
