# SPDX-License-Identifier: CECILL-2.1

module FiniteStateAutomata

using LinearAlgebra
using ChainRulesCore
using SparseArrays
using Semirings

#Base.oneunit(K::Type{<:Semiring}) = one(K)

export
    # concrete types
    FST,
    #DenseFST,

    # Accessors / properties
    α,
    T,
    ω,
    ρ,
    λ,
    #nstates,
    #nedges,
    #accessible,
    #coaccessible,

    # FST operations
    W,
    #Π₁,
    #Π₂,
    closure
    #determinize,
    #renorm
#    statemap,

#    connect,
#    determinize,
#    globalrenorm,
#    minimize,
#    propagate,
#    renorm

include("transitionmatrix.jl")
include("abstractfst.jl")
include("fst.jl")
#include("dense_fsa.jl")

include("ops.jl")

#include("totalweight.jl")
#include("project.jl")
#include("reverse.jl")
#include("union.jl")
#include("cat.jl")
#include("closure.jl")
#include("determinize.jl")
#include("renorm.jl")

#include("kron.jl")
#include("statemap.jl")

#include("intersect.jl")
#include("autograd.jl")

end

