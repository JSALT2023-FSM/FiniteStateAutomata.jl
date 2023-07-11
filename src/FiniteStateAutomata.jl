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

export vector2dict_lod, vector2dict_sod, dict2coo, dict2coo_csc

include("fst/fstconversion.jl")

#=====================================================================#
# Vizualisation.
#=====================================================================#

# export draw

# include("fst/graphviz.jl")

#=====================================================================#
# FST operations
#=====================================================================#

export dense_composition_sod, dense_composition_lod

include("fst/dense_ops.jl")

export sparse_composition_sod,  sparse_composition_lod
export sparse_coo_composition_lod
export sparse_vec_composition_lod_mt, sparse_coo_composition_lod_mt, sparse_coo2dict_composition_lod_mt

include("fst/sparse_ops.jl")

#=====================================================================#
# Concrete types.
#=====================================================================#

#=====================================================================#
# Loading/Saving FSTs.
#=====================================================================#

export draw, dot, loadsymbols, compile, Dot2SVG, Dot2PNG

include("fst/io.jl")


end

