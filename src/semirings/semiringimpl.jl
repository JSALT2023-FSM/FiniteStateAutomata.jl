# SPDX-License-Identifier: CECILL-2.1

# By default just return the the "val" property.
val(x::Semiring) = x.val

"""
    struct BoolSemiring <: Semiring
        val::Bool
    end

Boolean semiring: ``\\langle \\{0, 1\\}, \\lor, \\land, 0, 1 \\rangle``.
"""
struct BoolSemiring <: Semiring
    val::Bool
end

⊕(x::BoolSemiring, y::BoolSemiring) = BoolSemiring(x.val || y.val)
⊗(x::BoolSemiring, y::BoolSemiring) = BoolSemiring(x.val && y.val)
Base.zero(::Type{BoolSemiring}) = BoolSemiring(false)
Base.one(::Type{BoolSemiring}) = BoolSemiring(true)


"""
    struct LogSemiring{T,a} <: Semiring
        val::T
    end

Logarithmic semiring: ``\\langle \\mathbb{R} \\cup \\{- \\infty \\}, \\oplus_{\\log}, +, -\\infty, 0 \\rangle``
where

```math
x \\oplus y = \\frac{1}{a} \\log ( e^{ax} + e^{ay} ).
```
"""
struct LogSemiring{T<:AbstractFloat,a} <: Semiring
    val::T
end

⊕(x::LogSemiring{T,a}, y::LogSemiring{T,a}) where {T,a} = LogSemiring{T,a}(logaddexp(a*x.val, a*y.val)/a)
⊗(x::LogSemiring{T,a}, y::LogSemiring{T,a}) where {T,a} = LogSemiring{T,a}(x.val + y.val)
⊘(x::LogSemiring{T,a}, y::LogSemiring{T,a}) where {T,a} = LogSemiring{T,a}(x.val - y.val)
Base.inv(x::LogSemiring{T,a}) where {T,a} = LogSemiring{T,a}(-x.val)
Base.zero(K::Type{<:LogSemiring{T}}) where T = K(T(-Inf))
Base.one(K::Type{<:LogSemiring{T}}) where T = K(T(0))

"""
    struct ProbSemiring{T<:Real} <: Semiring{T}
        val::T
    end

Probability semiring ``\\langle (\\mathbb{R}_+``, +, \\cdot, 0, 1 \\rangle``.
"""
struct ProbSemiring{T<:AbstractFloat} <: Semiring
    val::T
end

⊕(x::ProbSemiring, y::ProbSemiring) = ProbSemiring(x.val + y.val)
⊗(x::ProbSemiring, y::ProbSemiring) = ProbSemiring(x.val * y.val)
⊘(x::ProbSemiring, y::ProbSemiring) = ProbSemiring(x.val / y.val)
Base.inv(x::ProbSemiring) = ProbSemiring(inv(x.val))
Base.zero(K::Type{<:ProbSemiring{T}}) where T = K(zero(T))
Base.one(K::Type{<:ProbSemiring{T}}) where T = K(one(T))

"""
    const TropicalSemiring{T} = LogSemiring{T,Inf} where T

Tropical semiring: ``\\langle \\mathbb{R} \\cup \\{- \\infty \\}, max, +, -\\infty, 0 \\rangle``.
"""
const TropicalSemiring{T} = LogSemiring{T,Inf} where T

⊕(x::TropicalSemiring{T}, y::TropicalSemiring{T}) where T = TropicalSemiring{T}(max(x.val, y.val))

