# SPDX-License-Identifier: CECILL-2.1

struct VectorFST{S,L} <: MutableFST{S,L}
    states::Vector{Vector{Tuple{Int, L, S}}}
    initstate::Int
    finalweights::Vector{S}
end


# AbstractFST operations
states(fst::VectorFST) = 1:length(fst.states)
arcs(fst::VectorFST, q) = fst.states[q]
initstate(fst::VectorFST) = fst.initstate
finalweight(fst::VectorFST, q) = fst.finalweights[q]


# MutableFST operations
setinitstate!(fst::VectorFST, q) = begin fst.initstate = q end

function setfinalstate!(fst::VectorFST{S}, q; w=one(S)) where S
    fst.finalweights[q] = w
end

function addstate!(fst::VectorFST{S}) where S
    push!(fst.states, [])
    push!(fst.finalweights, zero(S))
    length(fst.states)
end

function addarc!(fst::VectorFST, q, arc::Arc)
    push!(fst.states[q], arc)
    length(fst.states[q])
end

function deletestates!(fst::VectorFST, qs::AbstractVector)
    if fst.initstate ∈ qs
        error("Cannot remove initial state! Call `setinitstate` first.")
    end
    deleteat!(fst.states, qs)
    deleteat!(fst.finalweights, qs)
end

function deletearcs!(fst::VectorFST, q, arc_ids::AbstractVector)
    deleteat!(fst.states[q], arc_ids)
end

function deletearcs!(fst::VectorFST, q)
    fst.states[q] = []
end
