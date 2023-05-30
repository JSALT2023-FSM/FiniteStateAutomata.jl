# SPDX-License-Identifier: CECILL-2.1

#"""
#    accessible(A::AbstractFST{K}) where K
#
#Return a vector `x` of dimension `nstates(A)` where `x[i]` is `one(K)`
#if the state `i` is not accessible, `zero(K)` otherwise.
#"""
#function accessible(A::AbstractFST{K}) where K
#    vₙ = α(A)
#
#    m = ones(Bool, nstates(A))
#	m[findnz(vₙ)[1]] .= false
#
#	while nnz(vₙ) > 0
#		uₙ = T(A)' * vₙ
#		vₙ = uₙ .* m
#		m[findnz(uₙ)[1]] .= false
#	end
#	.! m
#end
#
#"""
#    coaccessible(A::AbstractFST{K}) where K
#
#Returns a vector `x` of dimension `nstates(A)` where `x[i]` is `one(K)`
#if the state `i` is not accessible, `zero(K)` otherwise.
#"""
#function coaccessible(A::AbstractFST)
#    vₙ = ω(A)
#
#    m = ones(Bool, nstates(A))
#	m[findnz(vₙ)[1]] .= false
#
#	while nnz(vₙ) > 0
#		uₙ = T(A) * vₙ
#		vₙ = uₙ .* m
#		m[findnz(uₙ)[1]] .= false
#	end
#	.! m
#end

"""
    W(A)

Return the total weight of `A`, i.e. the ``\\oplus``-sum of all the
path's weight in `A`.
"""
W(A::TransducerOrAcceptor) =
    ρ(A) ⊕ (transpose(α(A)) * MatrixPowerSum(T(A)) * ω(A)) # TODO optimize with dot(., ., .)

function closure(A::TransducerOrAcceptor; plus = false)
    K = eltype(α(A))
    TA = T(A)
    TB = TransitionMatrix(
        TA.S,
        hcat(TA.U, ω(A)),
        blockdiag(TA.E, sparse([1], [1], one(K), 1, 1)),
        vcat(TA.V, transpose(α(A)))
    )

    FST(α(A), TB, ω(A), ρ(A), λ(A))
end

#"""
#    cat(A1[, A2, ...])
#
#Return the concatenation of the given FSTs.
#"""
#function Base.cat(A::TransducerOrAcceptor, B::TransducerOrAcceptor)
#    K = promote_type(eltype(α(A)), eltype(α(B)))
#    TA, TB = T(A), T(B)
#    TC = TransitionMatrix(
#        TA.S,
#        hcat(TA.U, ω(A)),
#        blockdiag(TA.E, sparse([1], [1], one(K), 1, 1)),
#        vcat(TA.V, transpose(α(A)))
#    )
#end
#Base.cat(A1::AbstractFST{K}, AN::AbstractFST{K}...) where K = foldl(cat, AN, init = A1)

"""
    union(A1[, A2, ...])
    A1 ∪ A2

Return the union of the given FST.
"""
function Base.union(A::TransducerOrAcceptor, B::TransducerOrAcceptor)
    K = promote_type(eltype(α(A)), eltype(α(B)))
    TA, TB = T(A), T(B)
    TC = TransitionMatrix(
        blockdiag(TA.S, TB.S),
        blockdiag(TA.U, TB.U),
        blockdiag(TA.E, TB.E),
        blockdiag(TA.V, TB.V)
    )

    FST(vcat(α(A), α(B)), TC, vcat(ω(A), ω(B)), ρ(A) ⊕ ρ(B), vcat(λ(A), λ(B)))

end

Base.union(A1::AbstractFST{K}, AN::AbstractFST{K}...) where K = foldl(union, AN, init = A1)

#"""
#    connect(A::AbstractFST)
#
#Return a FST ``B`` equivalent of ``A`` such that all states in ``B``
#are accessible and coaccessible.
#"""
#function connect(A::AbstractFST{K}) where K
#    m = accessible(A) .* coaccessible(A)
#    I = findall(m)
#    M = sparse(I, 1:length(I), one(K), nstates(A), length(I))
#    FST(
#        M' * α(A),
#        M' * T(A) * M,
#        M' * ω(A),
#        ρ(A),
#        λ(A)[I]
#    )
#end
#
#@inline mergelabels(x) = x
#@inline mergelabels(x, y) = (x..., y...)
#@inline mergelabels(x, y, z...) = mergelabels(mergelabels(x, y), z...)
#
#function Base.replace(f::Function, M::AbstractFST, Ms)
#    R = ρ.(Ms)
#    I = findall(! iszero, R)
#    M = addskipedges(M, I, R[I])
#
#    A = blockdiag([α(m)[:,1:1] for m in Ms]...)
#    Ω = blockdiag([ω(m)[:,1:1] for m in Ms]...)
#    D = spdiagm([ρ(m) for m in Ms])
#    FST(
#        vcat([α(M)[i] * α(m) for (i, m) in enumerate(Ms)]...),
#        blockdiag([T(m) for m in Ms]...) + Ω * (T(M) + T(M) * D * T(M)') * A',
#        vcat([ω(M)[i] * ω(m) for (i, m) in enumerate(Ms)]...),
#        ρ(M),
#        vcat([f.([λ(M)[i]], λ(m)) for (i, m) in enumerate(Ms)]...)
#    )
#end
#
#function Base.replace(new::Function, M::AbstractFST, labelfn::Function)
#    replace(labelfn, M, [new(i) for i in 1:nstates(M)])
#end
#
#Base.replace(new::Function, M::AbstractFST) =
#    replace(new, M, mergelabels)
#
##= Functions for acyclic FST =#
#
#"""
#    propagate(A::AbstractAcyclicFST[; forward = true])
#
#Multiply the weight of each arc by the sum-product of all the prefix
#paths. If `forward` is `false` the arc's weight is multiply by the
#sum-product of all the suffix paths of this arc.
#"""
#function propagate(A::AbstractAcyclicFST; forward = true)
#    v = forward ? α(A) : ω(A)
#    M = forward ? T(A)' : T(A)
#    σ = v
#    while nnz(v) > 0
#        v = M * v
#        σ += v
#    end
#    σ
#
#    D = spdiagm(σ)
#    AcyclicFST(FST(σ .* α(A), T(A) * D, ω(A), ρ(A), λ(A)))
#end
#
#"""
#    globalrenorm(A::AbstractAcyclicFST)
#
#Global renormalization of the weights of `A`.
#"""
#globalrenorm(A::AbstractAcyclicFST) = propagate(A; forward = false) |> renorm
#
