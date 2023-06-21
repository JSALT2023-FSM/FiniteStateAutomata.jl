# SPDX identifier: CECILL-2.1

"""
    const ϵ::Int = 0

Epsilon label symbol represented as `0` by convention.
"""
const ϵ::Int = 0

"""
    const LabelSymbol = Integer

Label symbol.
"""
const LabelSymbol = Integer

"""
    const LabelSymbolMapping = Pair{<:LabelSymbol,<:LabelSymbol}

Label symbol fo an arc. It is stored as an integer.
"""
const LabelSymbolMapping = Pair{<:LabelSymbol,<:LabelSymbol}

"""
    const Label = Union{LabelSymbol,LabelSymbolMapping}

Arc's label. It can be either a unique symbol (for acceptor) or a
mapping rule `in_symbol => out_symbol` (for transducer).
"""
const Label = Union{LabelSymbol, LabelSymbolMapping}

isepsilon(l::LabelSymbol) = l == ϵ
isepsilon(l::LabelSymbolMapping) = isepsilon(first(l)) && isepsilon(last(l))

hasinputepsilon(l::LabelSymbol) = isepsilon(l)
hasoutputepsilon(l::LabelSymbol) = isepsilon(l)
hasinputepsilon(l::LabelSymbolMapping) = isepsilon(first(l))
hasoutputepsilon(l::LabelSymbolMapping) = isepsilon(last(l))

"""
    abstract type AbstractFST{S<:Semiring,L<:Label} end

Abstract base type for all FST. `S` is the weight semiring and `L` is
the label type. For acceptors, `L<:LabelSymbol` whereas for transducers
`L<:LabelMapping`.
"""
abstract type AbstractFST{S<:Semiring,L<:Label} end

"""
    M(fst)

Return a 3D tensor representing the arcs of the FST.
"""
M(::AbstractFST)


"""
    α(fst)

Return the vector of initial states of `A`.
"""
α(::AbstractFST)

"""
    ω(fst)

Return the vector of final states of `A`.
"""
ω(::AbstractFST)

"""
    semiring(fst)

Return the semiring type of `fst`.
"""
semiring(fst::AbstractFST{S}) where S = S

"""
    nstates(fst)

Return the number of states in `fst`.

See also: [`states`](@refs), [`arcs`](@ref), and [`narcs`](@ref).
"""
nstates(::AbstractFST)

"""
    narcs(fst)

Return the number of arcs in `fst`.

See also: [`arcs`](@ref), [`states`](@ref) and [`nstates`](@ref).
"""
narcs(::AbstractFST)

"""
    arcs(fst)

Iterator over the arcs of `fst`. Each element given by the iterator
is of the form `(src, dest, label, weight)`.

See also: [`narcs`](@ref), [`states`](@ref) and [`nstates`](@ref).
"""
arcs(fst::AbstractFST)

"""
    states(fst)

Iterator over the states of `fst`. Each element given by the iterator
is of the form `(stateid, initialweight, finalweight)`.

See also: [`nstates`](@ref), [`arcs`](@ref) and [`narcs`](@ref).
"""
states(fst::AbstractFST)

