# SPDX-License-Identifier: CECILL-2.1

module FiniteStateAutomata

include("semirings/Semirings.jl")
include("sparsesemimodules/SparseSemimodules.jl")

export SparseVector
export SparseMatrix

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

include("ops.jl")

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

