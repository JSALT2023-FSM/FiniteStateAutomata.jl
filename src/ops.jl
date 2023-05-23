# SPDX-License-Identifier: CECILL-2.1

"""
    accessible(A::AbstractFST{K}) where K

Return a vector `x` of dimension `nstates(A)` where `x[i]` is `one(K)`
if the state `i` is not accessible, `zero(K)` otherwise.
"""
function accessible(A::AbstractFST{K}) where K
    vₙ = α(A)

    m = ones(Bool, nstates(A))
	m[findnz(vₙ)[1]] .= false

	while nnz(vₙ) > 0
		uₙ = T(A)' * vₙ
		vₙ = uₙ .* m
		m[findnz(uₙ)[1]] .= false
	end
	.! m
end

"""
    coaccessible(A::AbstractFST{K}) where K

Returns a vector `x` of dimension `nstates(A)` where `x[i]` is `one(K)`
if the state `i` is not accessible, `zero(K)` otherwise.
"""
function coaccessible(A::AbstractFST)
    vₙ = ω(A)

    m = ones(Bool, nstates(A))
	m[findnz(vₙ)[1]] .= false

	while nnz(vₙ) > 0
		uₙ = T(A) * vₙ
		vₙ = uₙ .* m
		m[findnz(uₙ)[1]] .= false
	end
	.! m
end

"""
    connect(A::AbstractFST)

Return a FST ``B`` equivalent of ``A`` such that all states in ``B``
are accessible and coaccessible.
"""
function connect(A::AbstractFST{K}) where K
    m = accessible(A) .* coaccessible(A)
    I = findall(m)
    M = sparse(I, 1:length(I), one(K), nstates(A), length(I))
    FST(
        M' * α(A),
        M' * T(A) * M,
        M' * ω(A),
        ρ(A),
        λ(A)[I]
    )
end

@inline mergelabels(x) = x
@inline mergelabels(x, y) = (x..., y...)
@inline mergelabels(x, y, z...) = mergelabels(mergelabels(x, y), z...)

function Base.replace(f::Function, M::AbstractFST, Ms)
    R = ρ.(Ms)
    I = findall(! iszero, R)
    M = addskipedges(M, I, R[I])

    A = blockdiag([α(m)[:,1:1] for m in Ms]...)
    Ω = blockdiag([ω(m)[:,1:1] for m in Ms]...)
    D = spdiagm([ρ(m) for m in Ms])
    FST(
        vcat([α(M)[i] * α(m) for (i, m) in enumerate(Ms)]...),
        blockdiag([T(m) for m in Ms]...) + Ω * (T(M) + T(M) * D * T(M)') * A',
        vcat([ω(M)[i] * ω(m) for (i, m) in enumerate(Ms)]...),
        ρ(M),
        vcat([f.([λ(M)[i]], λ(m)) for (i, m) in enumerate(Ms)]...)
    )
end

function Base.replace(new::Function, M::AbstractFST, labelfn::Function)
    replace(labelfn, M, [new(i) for i in 1:nstates(M)])
end

Base.replace(new::Function, M::AbstractFST) =
    replace(new, M, mergelabels)

#= Functions for acyclic FST =#

"""
    propagate(A::AbstractAcyclicFST[; forward = true])

Multiply the weight of each arc by the sum-product of all the prefix
paths. If `forward` is `false` the arc's weight is multiply by the
sum-product of all the suffix paths of this arc.
"""
function propagate(A::AbstractAcyclicFST; forward = true)
    v = forward ? α(A) : ω(A)
    M = forward ? T(A)' : T(A)
    σ = v
    while nnz(v) > 0
        v = M * v
        σ += v
    end
    σ

    D = spdiagm(σ)
    AcyclicFST(FST(σ .* α(A), T(A) * D, ω(A), ρ(A), λ(A)))
end

"""
    globalrenorm(A::AbstractAcyclicFST)

Global renormalization of the weights of `A`.
"""
globalrenorm(A::AbstractAcyclicFST) = propagate(A; forward = false) |> renorm

