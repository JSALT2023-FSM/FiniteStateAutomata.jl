# SPDX-License-Identifier: CECILL-2.1

"""
    struct FST{S,L} <: AbstractFST{S,L}
        M::AbstractMatrix{S}
        α::AbstractVector{S}
        ω::AbstractVector{S}
    end

Generic Finite State Automaton.
"""
struct FST{S,L} <: AbstractFST{S,L}
    M::AbstractArray{S,3}
    α::AbstractVector{S}
    ω::AbstractVector{S}
    λ::AbstractVector{L}
end

M(fst::FST) = fst.M
α(fst::FST) = fst.α
ω(fst::FST) = fst.ω

