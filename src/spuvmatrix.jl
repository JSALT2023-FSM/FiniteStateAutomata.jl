# SPDX identifier: CECILL-2.1

struct MatrixPowerSum{K} <: AbstractMatrix{K}
    A::AbstractMatrix{K}
end

Base.parent(A::MatrixPowerSum) = A.A
Base.size(M::MatrixPowerSum) = size(M.S)

function Base.:*(xᵀ::Transpose{K, AbstractVector{K}}, M::MatrixPowerSum{K}) where K
    uₙᵀ = xᵀ
    i = 0
    while nnz(parent(uₙᵀ)) > 0
        i >= size(M, 1) || throw(ArgumentError("matrix `M` is not nilpotent"))
		uₙᵀ = uₙᵀ * T(A)
	end
    uₙᵀ
end

function Base.:*(M::MatrixPowerSum{K}, x::AbstractVector{K}) where K
    vₙ = x
    i = 0
    while nnz(vₙ) > 0
        i >= size(M, 1) || throw(ArgumentError("matrix `M` is not nilpotent"))
		vₙ = T(A) * vₙ
	end
    vₙ
end

#= Factorized matrix =#

struct FactorizedMatrix{K} <: AbstractMatrix{K}
    S::AbstractMatrix{K}
    U::AbstractMatrix{K}
    E::AbstractMatrix{K}
    V::AbstractMatrix{K}
end

Base.size(M::FactorizedMatrix) = size(M.S)

function Base.:*(xᵀ::Transpose{K, AbstractVector{K}}, M::FactorizedMatrix{K}) where K
    xᵀ * M.S + (xᵀ * M.U) * M.V
end

function Base.:*(M::FactorizedMatrix{K}, x::AbstractVector{K}) where K
    M.S * x + M.U * (M.V * x)
end

