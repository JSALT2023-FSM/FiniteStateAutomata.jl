# SPDX identifier: CECILL-2.1

const LabelSymbol = Union{<:Integer, NTuple{N,<:Integer} where N}
const LabelSymbolMapping = Pair{<:LabelSymbol,<:LabelSymbol}
const Label = Union{LabelSymbol, LabelSymbolMapping}

isepsilon(l::Integer) = l == 0
isepsilon(l::NTuple{1}) = isepsilon(l[1])
isepsilon(l::NTuple{N}) where N = isepsilon(l[1]) && isepsilon(l[2:end])
isepsilon(l::LabelSymbolMapping) = isepsilon(first(l)) && isepsilon(last(l))

isinputepsilon(l::LabelSymbol) = isepsilon(l)
isoutputepsilon(l::LabelSymbol) = isepsilon(l)
isinputepsilon(l::LabelSymbolMapping) = isepsilon(first(l))
isoutputepsilon(l::LabelSymbolMapping) = isepsilon(last(l))

"""
    abstract type AbstractFST{K<:Semiring,L<:Union{Label,LabelMapping}} end

Abstract base type for all FST. `K` is the weight
semiring and `L` is the label type or a label mapping. Acceptors is the
subset of FSTs for which `L` is a `Label` (or equivalently the label
mapping is just the identity function).
"""
abstract type AbstractFST{K<:Semiring,L<:Label} end

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
    semiring(fst)

Return the semiring type of `fst`.
"""
semiring(fst::AbstractFST{K}) where K = K

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

struct ArcIterator
    fst::AbstractFST
end

Base.length(it::ArcIterator) = narcs(it.fst)

function Base.iterate(it::ArcIterator, a = 1)
    a > narcs(it.fst)  && return nothing
    I, V1 = findnz(S(it.fst)[:, a])
    J, V2 = findnz(D(it.fst)[:, a])
    arc = (I[1], J[1], λ(it.fst)[a], V1[1] ⊗ V2[1])
    return (arc, a+1)
end

"""
    states(fst)

Iterator over the states of `fst`.
"""
arcs(fst::AbstractFST) = ArcIterator(fst)

struct StateIterator
    fst::AbstractFST
end

Base.length(it::StateIterator) = nstates(it.fst)

function Base.iterate(it::StateIterator, q = 1)
    q > nstates(it.fst)  && return nothing
    return ((q, α(it.fst)[q], ω(it.fst)[q]), q+1)
end

"""
    states(fst)

Iterator over the states of `fst`.
"""
states(fst::AbstractFST) = StateIterator(fst)

