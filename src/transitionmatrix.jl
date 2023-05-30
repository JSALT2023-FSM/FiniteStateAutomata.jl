# SPDX identifier: CECILL-2.1

struct MatrixPowerSum{K} <: AbstractMatrix{K}
    A::AbstractMatrix{K}
end

Base.parent(A::MatrixPowerSum) = M.A
Base.size(M::MatrixPowerSum) = size(M.A)

# TODO performance can be improved by preallocating buffer.
function Base.:*(xᵀ::Transpose{K,<:AbstractVector{K}}, M::MatrixPowerSum{K}) where K
    uₙᵀ = xᵀ
    acc = xᵀ
    n = 1
    while nnz(parent(uₙᵀ)) > 0
        n > size(M, 1) && throw(ArgumentError("matrix is not nilpotent"))
		uₙᵀ = uₙᵀ * M.A
        acc += uₙᵀ
        n += 1
	end
    acc
end

# TODO performance can be improved by preallocating buffer.
function Base.:*(M::MatrixPowerSum{K}, x::AbstractVector{K}) where K
    vₙ = x
    acc = x
    n = 1
    while nnz(vₙ) > 0
        n > size(M, 1) && throw(ArgumentError("matrix is not nilpotent"))
		vₙ = M.A * vₙ
        acc += vₙ
        n += 1
	end
    acc
end

#= Factorized matrix =#

struct TransitionMatrix{K} <: AbstractMatrix{K}
    S::AbstractMatrix{K}
    U::AbstractMatrix{K}
    E::AbstractMatrix{K}
    V::AbstractMatrix{K}
end

Base.size(M::TransitionMatrix) = size(M.S)
Base.:*(xᵀ::Transpose{K,<:AbstractVector}, M::TransitionMatrix{K}) where K =
    xᵀ * M.S + ((xᵀ * M.U) * M.E) * M.V

Base.:*(M::TransitionMatrix{K}, x::AbstractVector{K}) where K =
    M.S * x + M.U * (M.E * (M.V * x))

