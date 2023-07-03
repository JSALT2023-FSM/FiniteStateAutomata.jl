# SPDX identifier: CECILL-2.1

"""
    const Arc{S} = Tuple{State, Label, Label, S} where {S}

Arc type, store the destination state, the input and output label and
a weight.
"""
const Arc{S} = Tuple{State, Label, Label, S} where {S}

#=====================================================================#
# Abstact FST interface.
#=====================================================================#

"""
    abstract type AbstractFST{S<:Semiring,L<:Label} end

Abstract base type for all FST. `S` is the weight semiring and `L` is
the label type.
"""
abstract type AbstractFST{S<:Semiring} end

"""
    semiring(fst)

Return the semiring type of `fst`.
"""
semiring(fst::AbstractFST{S}) where S = S

"""
    numarcs(fst, q)

Return the number of arcs leaving state `q` in `fst`.
"""
numarcs(fst::AbstractFST, q) = length(arcs(fst, q))

"""
    states(fst)

Iterator over the states of `fst`.
"""
states(fst::AbstractFST)

"""
    arcs(fst, q)

Iterator over the arcs leaving state `q` in `fst`.
"""
arcs(fst::AbstractFST, q)

"""
    initstate(fst)

Return the initial state.
"""
initstate(fst::AbstractFST)

"""
    isinit(fst, q)

Check whether state `q` is initial.
"""
isinit(fst::AbstractFST, q) = q == initstate(fst)

"""
    finalweight(fst, q)

Return the weight of a final state if `q` is final, 0̄ otherwise.
"""
finalweight(fst::AbstractFST, q)

"""
    isfinal(fst, q)

Returns true if `q` is final in `fst`.
"""
isfinal(fst::AbstractFST, q) = !iszero(finalweight(fst, q))

"""
    finalstates(fst)

Returns iterator over final states.
"""
finalstates(fst::AbstractFST) = filter(q -> isfinal(fst, q), states(fst))

#=====================================================================#
# Expanded FST interface.
#=====================================================================#

"""
    abstract type ExpandedFST{S} <: AbstractFST{S} end

Base class for FST expanded in memory.
"""
abstract type ExpandedFST{S} <: AbstractFST{S} end

"""
    numstates(fst)

Return the number of states in `fst`.
"""
numstates(fst::ExpandedFST) = length(states(fst))

#=====================================================================#
# Mutable FST interface.
#=====================================================================#

"""
    abstract type MutableFST{S} <: ExpandedFST{S} end

Base class for modifiable FST.
"""
abstract type MutableFST{S} <: ExpandedFST{S} end

"""
    addstate!(fst, ...)

Add empty state to `fst`
Returns state ID.
"""
addstate!(fst::MutableFST)

"""
    setinitstate!(fst, q)

Set state `q` as initial state in `fst`.
"""
setinitstate!(fst::MutableFST, q)

"""
    setfinalstate!(fst, q, w)

Set state `q` as final state with weight `w`.
"""
setfinalstate!(fst::MutableFST{S}, q, weight::S) where S

"""
    addarc!(fst, q, arc)

Add tranisiton `arc` leading from state `q`.
Returns arc ID.
"""
addarc!(fst::MutableFST{S}, q, a::Arc{S}) where {S}

"""
    deletestate!(fst, q)

Delete state `q` from `fst`.
"""
deletestate!(fst::MutableFST, q)

"""
    deletestates!(fst, states)

Delete states `qs` from `fst`.
"""
function deletestates!(fst::MutableFST, qs)
    for state in qs
        deletestate!(fst, q)
    end
    fst
end

"""
    deletestates!(fst)

Delete all states from `fst`.
"""
deletestates!(fst::MutableFST) = deletestates!(fst, states(fst))

"""
    deletearcs!(fst, q, arc_ids)

Delete arcs leading from state `q`.
"""
deletearcs!(fst::MutableFST, q, arc_ids)

"""
    deletearc!(fst, q, arc_id)

Delete arc leading from state `q`.
"""
deletearc!(fst::MutableFST, q, arc_id) = deletearcs!(fst, q, [arc_id])

"""
    deletearcs!(fst, q)

Delete all arcs leading from the state `q`.
"""
deletearcs!(fst::MutableFST, q)

