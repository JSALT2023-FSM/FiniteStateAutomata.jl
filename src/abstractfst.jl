# SPDX identifier: CECILL-2.1

const Label = Int
const LabelMapping = Pair{<:Label,<:Label}

"""
    abstract type AbstractFST{K<:Semiring,L<:Union{Label,LabelMapping}} end

Abstract base type for all FST. `K` is the weight
semiring and `L` is the label type or a label mapping. Acceptors is the
subset of FSTs for which `L` is a `Label` (or equivalently the label
mapping is just the identity function).
"""
abstract type AbstractFST{K<:Semiring,L<:Union{Label,LabelMapping}} end

const Acceptor = AbstractFST{<:Semiring,<:Label}

"""
    α(A)

Return the vector of initial states of `A`.
"""
α(::AbstractFST)

"""
    S(A)

Return the source incidence matrix of `A`.
"""
S(::AbstractFST)

"""
    D(A)

Return the destination incidence matrix of `A`.
"""
D(::AbstractFST)

"""
    ω(A)

Return the vector of final states of `A`.
"""
ω(::AbstractFST)

"""
    λ(A)

Return the states' label of `A`.
"""
λ(::AbstractFST)

"""
    nstates(A)

Return the number of states in `A`.
"""
nstates(A::AbstractFST) = size(S(A), 1)

"""
    narcs(A)

Return the number of arcs in `A`.
"""
narcs(A::AbstractFST) = size(S(A), 2)


"""
    isymbols(A)

Return the input symbol table of `A`.
"""
isymbols(A::AbstractFST)

"""
    osymbols(A)

Return the output symbol table of `A`.
"""
osymbols(A::AbstractFST)

