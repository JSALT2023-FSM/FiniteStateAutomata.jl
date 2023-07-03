# SPDX-License-Identifier: CECILL-2.1

Base.convert(::Type{TensorFST{S, T}}, vfst::VectorFST{S}) where {S, T<:Array{S, 4}} = densefst(vfst)
