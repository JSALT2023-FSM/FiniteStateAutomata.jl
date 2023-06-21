# SPDX-License-Identifier: CECILL-2.1

module FiniteStateAutomata

using Semirings

#=====================================================================#
# Abstract FST types and generic properties.
#=====================================================================#

export
    M,
    α,
    ω,
    arcs,
    narcs,
    nstates,
    states,
    semiring


include("abstractfst.jl")

#=====================================================================#
# Vizualisation.
#=====================================================================#

#include("graphviz.jl")
#include("io.jl")

#=====================================================================#
# FST operations
#=====================================================================#

#include("ops.jl")

#=====================================================================#
# Concrete types.
#=====================================================================#

export FST

include("fst.jl")

end

