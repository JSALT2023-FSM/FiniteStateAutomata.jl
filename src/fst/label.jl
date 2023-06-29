# SPDX identifier: CECILL-2.1

"""
    const SymbolId = Integer

Label's symbol identifier.
"""
const SymbolId = Integer

"""
    const SymbolId = Pair{<:SymbolId,<:SymbolId}

Symbol id for an arc.
"""
const Mapping = Pair{<:SymbolId,<:SymbolId}

"""
    const Label = Union{SymbolId,SymbolMapping}

Arc's label. It can be either a unique symbol for acceptor or a
mapping rule `in_symbol => out_symbol` for transducer.
"""
const Label = Union{SymbolId, Mapping}

isepsilon(l::SymbolId) = l == ϵ
isepsilon(l::Mapping) = isepsilon(first(l)) && isepsilon(last(l))

hasinputepsilon(l::SymbolId) = isepsilon(l)
hasinputepsilon(l::Mapping) = isepsilon(first(l))

hasoutputepsilon(l::SymbolId) = isepsilon(l)
hasoutputepsilon(l::Mapping) = isepsilon(last(l))

inputsymbol(l::SymbolId) = l
inputsymbol(l::Mapping) = first(l)
outputsymbol(l::SymbolId) = l
outputsymbol(l::Mapping) = last(l)
