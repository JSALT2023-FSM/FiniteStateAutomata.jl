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
    AbstractFST,
    Acceptor,

    # Concrete types
    Label,
    FST,

    # Accessors / properties
    α,
    S,
    D,
    ω,
    ρ,
    λ,
    nstates,
    narcs,
    arcs,
    states,
    semiring,
    accessible,
    coaccessible,

    powerseries,

    # Utilities
    compile,
    symboltable,
    print,
    draw,

    # Unary FST operations
    filterarcs,
    filterstates,
    project,
    relabel,
    connect,
    closure,

    # FST operations
    W
    #determinize,
    #renorm
    #statemap,



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
include("io.jl")


end

