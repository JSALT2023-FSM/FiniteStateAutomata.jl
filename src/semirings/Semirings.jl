# SPDX-License-Identifier: CECILL-2.1

# Necessary for the logarithmic semiring
import LogExpFunctions: logaddexp, log1pexp

export
    # Types
    Semiring,
    BoolSemiring,
    LogSemiring,
    ProbSemiring,
    TropicalSemiring,
    ProductSemiring,

    # Semiring operations
    ⊕,
    ⊗,
    ⊘,

    # Get the "natural" value of a semiring element
    val

include("semiringtype.jl")
include("semiringimpl.jl")

# Allow semiring numbers to be used with common notation.
Base.:+(x::Semiring, y::Semiring) = x ⊕ y
Base.:*(x::Semiring, y::Semiring) = x ⊗ y
Base.:/(x::Semiring, y::Semiring) = x ⊘ y

# Real number operation intepreted as semiring operations
⊕(x::Real, y::Real) = x + y
⊗(x::Real, y::Real) = x * y
