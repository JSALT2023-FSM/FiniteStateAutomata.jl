# SPDX identifier: CECILL-2.1

struct MatrixPowerSeries{S} <: AbstractMatrix{S}
    A::AbstractMatrix{S}
end

powerseries(A::AbstractMatrix) = MatrixPowerSeries(A)

Base.parent(A::MatrixPowerSeries) = M.A
Base.size(M::MatrixPowerSeries) = size(M.A)

function Base.:*(xᵀ::Transpose{K,<:AbstractVector{K}}, M::MatrixPowerSeries{K}) where K
    # Vector of distances
    d = xᵀ

    uₙᵀ = xᵀ

    n = 1
    while nnz(parent(uₙᵀ)) > 0
        @show n
        n > size(M, 1) && throw(ArgumentError("matrix is not nilpotent"))
        uₙᵀ = uₙᵀ * M.A
        #d += uₙᵀ
        n += 1
        @show nnz(parent(uₙᵀ))
        println(uₙᵀ)
	end
    d
end

