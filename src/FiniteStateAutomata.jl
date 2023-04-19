# SPDX-License-Identifier: CECILL-2.1

module FiniteStateAutomata

using LinearAlgebra
using ChainRulesCore
using SparseArrays
using Semirings

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
    statemap,
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

include("kron.jl")
include("statemap.jl")
include("reverse.jl")
include("sum.jl")

include("ops.jl")

include("intersect.jl")
include("autograd.jl")

end

