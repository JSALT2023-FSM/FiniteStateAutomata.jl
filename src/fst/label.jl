# SPDX identifier: CECILL-2.1

"""
    const Label = Int

Label's symbol identifier.
"""
const Label = Int

isepsilon(l::Label) = l == ϵ

inputsymbol(l::SymbolId) = l
inputsymbol(l::Mapping) = first(l)
outputsymbol(l::SymbolId) = l
outputsymbol(l::Mapping) = last(l)
