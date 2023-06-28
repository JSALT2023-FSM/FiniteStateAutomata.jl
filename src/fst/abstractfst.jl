# SPDX identifier: CECILL-2.1

"""
    abstract type AbstractFST{S<:Semiring,L<:Label} end

Abstract base type for all FST. `S` is the weight semiring and `L` is
the label type.
"""
abstract type AbstractFST{S<:Semiring,L<:Label} end

"""
    semiring(fst)

Return the semiring type of `fst`.
"""
semiring(fst::AbstractFST{S}) where S = S

"""
    numstates(fst)

Return the number of states in `fst`.
"""
numstates(fst::AbstractFST) = length(states(fst))

"""
    numarcs(fst, q)

Return the number of arcs leaving state `q` in `fst`.
"""
numarcs(fst::AbstractFST, q) = length(arcs(fst, q))

"""
    states(fst)

Iterator over the states of `fst`.
"""
states(::AbstractFST)

"""
    arcs(fst, q)

Iterator over the arcs leaving state `q` in `fst`.
"""
arcs(fst::AbstractFST, q)

"""
    initstate(fst)

Return the initial state.
"""
initstate(::AbstractFST)

# M:
# Start -> s
# Final(state) -> w
# StateIterator ->
# Ar[
# -1 -> not label
