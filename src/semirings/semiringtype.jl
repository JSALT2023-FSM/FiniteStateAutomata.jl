# SPDX-License-Identifier: CECILL-2.1

"""
    abstract type Semiring <: Number end


Abstract type for a semiring ``\\langle S, \\oplus, \\otimes, \\bar{0}, \\bar{1} \\rangle``.

Subtypes of ``Semiring`` should implement the following functions:
- [`⊕`](@ref)
- [`⊗`](@ref)
- [`⊘`](@ref) (for divisible semiring only)
- [`Base.one]`
- [`Base.zero]`
"""
abstract type Semiring <: Number end

"""
    ⊕(x::Semiring, y::Semiring)
    x ⊕ y

Semiring addition.
"""
⊕

"""
    ⊗(x::Semiring, y::Semiring)
    x ⊗ y

Semiring multiplication.
"""
⊗

"""
    val(x::Semiring)

Return the "real" value / object wrapped in the semiring type
"""
val

# Patch Julia's "zero bug".
Base.zero(x::T) where T<:Semiring = zero(T)

function Base.show(io::IO, x::Semiring)
    if get(io, :compact, false) || get(io, :limit, false)
        print(io, val(x))
    else
        print(io, typeof(x), "(", val(x), ")")
    end
    io
end

