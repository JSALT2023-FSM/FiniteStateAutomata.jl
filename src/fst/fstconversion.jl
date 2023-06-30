# SPDX-License-Identifier: CECILL-2.1

Base.convert(::Type{TensorFST{S}}, vfst::VectorFST{S}) where S = densefst(vfst)
