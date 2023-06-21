# SPDX-License-Identifier: CECILL-2.1

module FiniteStateAutomata

using Semirings
using SparseSemimodules

#=====================================================================#
# Abstract FST types and generic properties.
#=====================================================================#

export
    M,
    α,
    ω,
    λ,
    arcs,
    narcs,
    nstates,
    states,
    semiring


include("abstractfst.jl")

#=====================================================================#
# Vizualisation.
#=====================================================================#

export draw

include("graphviz.jl")

#=====================================================================#
# FST operations
#=====================================================================#

#include("ops.jl")

#=====================================================================#
# Concrete types.
#=====================================================================#

export SparseFST

include("fst.jl")

#=====================================================================#
# Loading/Saving FSTs.
#=====================================================================#

export compile

include("io.jl")

end

