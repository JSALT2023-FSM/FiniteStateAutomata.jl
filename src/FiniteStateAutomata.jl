# SPDX-License-Identifier: CECILL-2.1

module FiniteStateAutomata

using LinearAlgebra
using SparseArrays

export
    # Abstract types
    AbstractFSA,
    AbstractAcyclicFSA,

    # concrete types
    AcyclicFSA,
    FSA,
    DenseFSA,
    IntersectedFSA,
    IntersectedDenseFSA,

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
    connect,
    determinize,
    globalrenorm,
    minimize,
    propagate,
    renorm

include("abstractfsa.jl")
include("fsa.jl")
include("ops.jl")
include("intersect.jl")
include("dense_fsa.jl")

end
