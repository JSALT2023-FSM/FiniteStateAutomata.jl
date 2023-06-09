# SPDX identifier: CECILL-2.1

struct MatrixPowerSeries{K} <: AbstractMatrix{K}
    A::AbstractMatrix{K}
    fn::Function
end

powerseries(A::AbstractMatrix) = powerseries(identity, A)
powerseries(fn::Function, A::AbstractMatrix) = MatrixPowerSeries(A, fn)

Base.parent(A::MatrixPowerSeries) = M.A
Base.size(M::MatrixPowerSeries) = size(M.A)

# TODO performance can be improved by preallocating buffer.
function Base.:*(xᵀ::Transpose{K,<:AbstractVector{K}}, M::MatrixPowerSeries{K}) where K
    uₙᵀ = M.fn(xᵀ)
    acc = xᵀ
    n = 1
    while nnz(parent(uₙᵀ)) > 0
        n > size(M, 1) && throw(ArgumentError("matrix is not nilpotent"))
        uₙᵀ = M.fn(uₙᵀ * M.A)
        acc += uₙᵀ
        n += 1
	end
    acc
end

# TODO performance can be improved by preallocating buffer.
function Base.:*(M::MatrixPowerSeries{K}, x::AbstractVector{K}) where K
    vₙ = M.fn(x)
    acc = x
    n = 1
    while nnz(vₙ) > 0
        n > size(M, 1) && throw(ArgumentError("matrix is not nilpotent"))
        vₙ = M.fn(M.A * vₙ)
        acc += vₙ
        n += 1
	end
    acc
end

