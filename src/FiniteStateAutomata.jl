# SPDX-License-Identifier: CECILL-2.1

module FiniteStateAutomata

using LinearAlgebra
using ChainRulesCore

# Currently we are using the standard Julia library of SparseArrays
# which supports only CPU. We will move to "SparseSemimodules` soon.
using SparseArrays

using Semirings

export
    # Abstract types
    Acceptor,

    # Concrete types
    FST,

    # Accessors / properties
    α,
    T,
    ω,
    ρ,
    λ,
    nstates,
    narcs,
    isymbols,
    osymbols,
    #accessible,
    #coaccessible,

    # Utilities
    compile,
    symboltable,
    print,

    # FST operations
    W,
    Π₁,
    Π₂,
    closure,
    #determinize,
    #renorm
    #statemap,

    draw

#    connect,
#    determinize,
#    globalrenorm,
#    minimize,
#    propagate,

# TODO: the TransitionMatrix structure and related types should
# move to `SparseSemimodules`.
include("transitionmatrix.jl")

include("abstractfst.jl")
include("graphviz.jl")
include("fst.jl")
include("ops.jl")

include("utils.jl")


end

