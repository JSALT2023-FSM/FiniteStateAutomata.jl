# SPDX-License-Identifier: CECILL-2.1

module FiniteStateAutomata

using LinearAlgebra
using SparseArrays

export
    # Abstract types
    AbstractFST,
    AbstractAcyclicFST,

    # concrete types
    AcyclicFST,
    FST,
    DenseFST,
    IntersectedFST,
    IntersectedDenseFST,

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

    # FST operations
    addskipedges,
    closure,
    connect,
    determinize,
    globalrenorm,
    minimize,
    propagate,
    renorm

include("abstractfst.jl")
include("iterate.jl")

include("fst.jl")
include("dense_fsa.jl")

include("reverse.jl")
include("ops.jl")
include("intersect.jl")

end

