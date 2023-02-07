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
    accessible,
    coaccessible,
    nstates,
    nedges,
    initstates,
    edges,
    finalstates,

    # FSA operations
    addskipedges,
    closure,
    globalrenorm,
    propagate,
    renorm

include("abstractfsa.jl")
include("fsa.jl")
include("ops.jl")

end
