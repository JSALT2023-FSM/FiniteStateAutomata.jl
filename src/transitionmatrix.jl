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
        i >= size(M, 1) || throw(ArgumentError("matrix is not nilpotent"))
		uₙᵀ = uₙᵀ * T(A)
	end
    uₙᵀ
end

function Base.:*(M::MatrixPowerSum{K}, x::AbstractVector{K}) where K
    vₙ = x
    i = 0
    while nnz(vₙ) > 0
        i >= size(M, 1) || throw(ArgumentError("matrix is not nilpotent"))
		vₙ = T(A) * vₙ
	end
    vₙ
end

#= Factorized matrix =#

# Returns a `Q` x `P` sparse matrix from a list of arcs.
# An arc is defined by a triplet `(i, j, v)`. `K` is element type
# of the matrix.
function _spm_from_list(K, Q, P, arcs)
    I, J, V = Int[], Int[], K[]
    for (src, dest, weight) in arcs
        push!(I, src)
        push!(J, dest)
        push!(V, weight)
    end
    sparse(I, J, V, Q, P)
end

struct TransitionMatrix{K} <: AbstractMatrix{K}
    S::AbstractMatrix{K}
    U::AbstractMatrix{K}
    E::AbstractMatrix{K}
    V::AbstractMatrix{K}
end

function TransitionMatrix(K, Q, direct_arcs, in_factors, factors, out_factors)
    S = _spm_from_list(K, Q, Q, direct_arcs)

    # Number of factors.
    P = max(
        maximum(t -> t[2], in_factors),
        maximum(t -> max(t[1], t[2]), factors),
        maximum(t -> t[1], out_factors)
    )

    U = _spm_from_list(K, Q, P, in_factors)
    E = _spm_from_list(K, P, P, factors)
    V = _spm_from_list(K, P, Q, out_factors)

    TransitionMatrix(S, U, E, V)
end

Base.size(M::TransitionMatrix) = size(M.S)

Base.:*(xᵀ::Transpose{K, AbstractVector{K}}, M::TransitionMatrix{K}) where K =
    xᵀ * M.S + (xᵀ * M.U) * M.V

Base.:*(M::TransitionMatrix{K}, x::AbstractVector{K}) where K =
    M.S * x + M.U * (M.V * x)

