# SPDX-License-Identifier: CECILL-2.1

module FiniteStateAutomata

using FileIO

#=====================================================================#
# Semirings definition.
#=====================================================================#

include("semirings/Semirings.jl")

#=====================================================================#
# Sparse linear algebra operations.
#=====================================================================#

include("sparsesemimodules/SparseSemimodules.jl")

export SparseMatrixCOO
# export SparseVector
# export SparseMatrix

#=====================================================================#
# Constants.
#=====================================================================#

# Does not export anything but sets the behavior.
include("fst/constants.jl")

#=====================================================================#
# Generic FST interface.
#=====================================================================#

export Arc
export arcs, numarcs, numstates, states, semiring
export initstate, isinit, finalweight, isfinal, finalstates
export addstate!, addarc!, setinitstate!, setfinalstate!
export deletestates!, deletestate!, deletearcs!, deletearc!

include("fst/abstractfst.jl")

#=====================================================================#
# Concrete FST implementations
#=====================================================================#

export VectorFST
export TensorFST

export M, α, ω

include("fst/vectorfst.jl")
include("fst/tensorfst.jl")

export vector2dict_lfo, vector2dict_sfo, dict2coo

include("fst/fstconversion.jl")

#=====================================================================#
# Vizualisation.
#=====================================================================#

# export draw

# include("fst/graphviz.jl")

#=====================================================================#
# FST operations
#=====================================================================#

export dense_composition_sfo, dense_composition_lfo, sparse_composition_lfo, sparse_composition_lfo

include("fst/dense_ops.jl")

export sparse_composition_sfo

include("fst/sparse_ops.jl")

#=====================================================================#
# Concrete types.
#=====================================================================#

#=====================================================================#
# Loading/Saving FSTs.
#=====================================================================#

export draw, dot, loadsymbols, compile, Dot2SVG

include("fst/io.jl")


end

