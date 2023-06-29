# SPDX-License-Identifier: CECILL-2.1

"""
    mutable struct VectorFST{S,L} <: MutableFST{S,L}
        states::Vector{Vector{Tuple{Int, L, S}}}
        initstate::Int
        finalweights::Vector{S}
    end

General-purpose mutable FST similar to the [*VectorFst*](https://www.openfst.org/twiki/bin/view/FST/FstAdvancedUsage#FST%20Types)
in Openfst.
"""
mutable struct VectorFST{S,L} <: MutableFST{S,L}
    states::Vector{Vector{Tuple{Int, L, S}}}
    initstate::Int
    finalweights::Vector{S}

    function VectorFST(states::AbstractVector{<:AbstractVector{Tuple{Int,L,S}}},
                            initstate,
                            finalweights) where {S,L}
        length(finalweights) == length(states) || throw(ArgumentError())
        new{S,L}(states, initstate, finalweights)
    end
end

# AbstractFST operations
states(fst::VectorFST) = 1:length(fst.states)
arcs(fst::VectorFST, q) = fst.states[q]
initstate(fst::VectorFST) = fst.initstate
finalweight(fst::VectorFST, q) = fst.finalweights[q]

# MutableFST operations
setinitstate!(fst::VectorFST, q) = fst.initstate = q

setfinalstate!(fst::VectorFST{S}, q; w=one(S)) where S = fst.finalweights[q] = w

function addstate!(fst::VectorFST{S}) where S
    push!(fst.states, [])
    push!(fst.finalweights, zero(S))
    length(fst.states)
end

function addarc!(fst::VectorFST, q, arc)
    push!(fst.states[q], arc)
    length(fst.states[q])
end

function deletestate!(fst::VectorFST, q)
    if fst.initstate == q
        error("Cannot remove initial state! Call `setinitstate!` first.")
    end
    deleteat!(fst.states, q)
    deleteat!(fst.finalweights, q)

    for qn in 1:length(fst.states)
        arcs = filter!(a -> a[1] != q, fst.states[qn])
        fst.states[qn] = [(d > q ? d - 1 : d, l, w) for (d, l, w) in arcs]
    end
end

deletearcs!(fst::VectorFST, q, arc_ids::AbstractVector) = deleteat!(fst.states[q], arc_ids)
deletearcs!(fst::VectorFST, q) = fst.states[q] = []

