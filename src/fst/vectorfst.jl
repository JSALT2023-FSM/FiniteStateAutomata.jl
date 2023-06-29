# SPDX-License-Identifier: CECILL-2.1

"""
    mutable struct VectorFST{S} <: MutableFST{S}
        arcs::Vector{Vector{Arc{S}}}
        initstate::Int
        finalweights::Vector{S}
    end

General-purpose mutable FST similar to the [*VectorFst*](https://www.openfst.org/twiki/bin/view/FST/FstAdvancedUsage#FST%20Types)
in Openfst.
"""
mutable struct VectorFST{S} <: MutableFST{S}
    arcs::AbstractVector{<:AbstractVector{Arc{S}}} 
    initstate::Int
    finalweights::AbstractVector{S}

    function VectorFST(arcs::AbstractVector{<:AbstractVector{Arc{S}}},
                            initstate,
                            finalweights) where {S}
        length(finalweights) == length(arcs) || throw(ArgumentError())
        new{S}(arcs, initstate, finalweights)
    end
end

# AbstractFST operations
states(fst::VectorFST) = 1:length(fst.arcs)
arcs(fst::VectorFST, q) = fst.arcs[q]
initstate(fst::VectorFST) = fst.initstate
finalweight(fst::VectorFST, q) = fst.finalweights[q]

# MutableFST operations
setinitstate!(fst::VectorFST, q) = fst.initstate = q

setfinalstate!(fst::VectorFST{S}, q; w=one(S)) where S = fst.finalweights[q] = w

function addstate!(fst::VectorFST{S}) where S
    push!(fst.arcs, [])
    push!(fst.finalweights, zero(S))
    length(fst.arcs)
end

function addarc!(fst::VectorFST, q, arc)
    push!(fst.arcs[q], arc)
    length(fst.arcs[q])
end

function deletestate!(fst::VectorFST, q)
    if fst.initstate == q
        error("Cannot remove initial state! Call `setinitstate!` first.")
    end
    deleteat!(fst.arcs, q)
    deleteat!(fst.finalweights, q)

    for qn in 1:length(fst.arcs)
        arcs = filter!(a -> a[1] != q, fst.arcs[qn])
        fst.arcs[qn] = [(d > q ? d - 1 : d, il, ol, w) for (d, il, ol, w) in arcs]
    end
end

deletearcs!(fst::VectorFST, q, arc_ids::AbstractVector) = deleteat!(fst.arcs[q], arc_ids)
deletearcs!(fst::VectorFST, q) = fst.arcs[q] = []

