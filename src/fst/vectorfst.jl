# SPDX-License-Identifier: CECILL-2.1

struct VectorFST{S,L} <: AbstractFST{S,L}
    states::Vector{Vector{L,S}}
end
