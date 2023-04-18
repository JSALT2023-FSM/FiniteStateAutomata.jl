# SPDX-License-Identifier: CECILL-2.1

function Base.intersect(A::AbstractFST, B::AbstractFST)
    kC = kron(A, B) do lx, ly
        lx == ly ? lx : "ϵ"
    end

    filter(kC) do q
        λ(kC)[q] != "ϵ"
    end
end

# Generic intersection
#struct IntersectedFST{K, L, T_A<:AbstractFST{K,L}, T_B<:AbstractFST{K,L}} <: AbstractFST{K, L}
#    A::T_A
#    B::T_B
#    M::AbstractMatrix{K}
#    λ::AbstractVector{L}
#end
#
#_computeM(::Type{K}, λ₁::AbstractVector{L}, λ₂::AbstractVector{L}) where {K,L} = begin
#    rows = []
#    labels = L[]
#    q = 0
#    for (i, s1) in enumerate(λ₁)
#        for (j, s2) in enumerate(λ₂)
#            if s1 == s2
#                q += 1
#                push!(rows, (i-1) * length(λ₂) + j)
#                push!(labels, s1)
#            end
#        end
#    end
#    sparse(rows, 1:q, one(K), length(λ₁) * length(λ₂), q), labels
#end
#
#IntersectedFST(A::AbstractFST{K}, B::AbstractFST{K}) where K = begin
#    M, labels = _computeM(K, λ(A), λ(B))
#    IntersectedFST(A, B, M, labels)
#end
#
#Base.intersect(A::AbstractFST, B::AbstractFST) = IntersectedFST(A, B)
#
#α(I::IntersectedFST) = I.M' * kron(α(I.A), α(I.B))
#T(I::IntersectedFST) = I.M' * kron(T(I.A), T(I.B)) * I.M
#ω(I::IntersectedFST) = I.M' * kron(ω(I.A), ω(I.B))
#ρ(I::IntersectedFST) = ρ(I.A) * ρ(I.B)

# Intersection with DenseFST
struct IntersectedDenseFST{K, L, T<:AbstractFST{K, L}} <: AbstractAcyclicFST{K, L}
    A::DenseFST{K, L}
    B::T
    C::AbstractMatrix{K}  # of shape |Σ| x nstates(A)
end

_computeC(::Type{K}, Σ::AbstractVector, labels::AbstractVector) where K = begin
    rows = []
    cols = []
    for (i,s) in enumerate(Σ)
        for (j, l) in enumerate(labels)
            if s == l
                push!(rows, i)
                push!(cols, j)
            end
        end
    end
    sparse(rows, cols, one(K), length(Σ), length(labels))
end

IntersectedDenseFST(A::DenseFST{K}, B::AbstractFST{K}) where K = begin
    C = _computeC(K, A.Σ, λ(B)) # TODO switch cols and rows
    IntersectedDenseFST(A, B, C)
end

Base.intersect(A::DenseFST, B::AbstractFST) = IntersectedDenseFST(A, B)
Base.intersect(A::AbstractFST, B::DenseFST) = IntersectedDenseFST(B, A) # we assume ∩ to be commutative w.r.t any K

α(I::IntersectedDenseFST) = begin
    h1 = I.A.H[:, 1]
    Q = nstates(I.B) * (size(I.A.H, 2) - 1)
    vals = (I.C' * h1) .* α(I.B)
    sparsevec(findnz(vals)..., Q)
end

ω(I::IntersectedDenseFST) = begin
    h = I.A.H[:, end]
    Q = nstates(I.B) * (size(I.A.H, 2) - 1)
    vals = (I.C' * h) .* ω(I.B)
    I, V = findnz(vals)
    sparsevec(I .+ (Q - length(vals)), V, Q)
end

T(I::IntersectedDenseFST) = begin
    H = I.A.H
    Tb = T(I.B)
    C = I.C
    K = eltype(H)
    Σ = length(I.A.Σ)
    Λ = nstates(I.B)
    Q = Λ * (size(H, 2) - 1)

    rows = []
    cols = []
    vals = K[]

    for n in 2:size(H, 2) - 1
        tmp = ones(K, Λ) * (C' * H[:, n])' .* Tb
        r, c, v = findnz(tmp)
        append!(rows, (n - 2) * Λ .+ r)
        append!(cols, (n - 1) * Λ .+ c)
        append!(vals, v)
    end
    sparse(rows, cols, vals, Q, Q)
end

λ(I::IntersectedDenseFST) = repeat(λ(I.B), size(I.A.H, 2) - 1)
ρ(I::IntersectedDenseFST) = ρ(I.B) * ρ(I.A)

function forward_backward(I::IntersectedDenseFST)
    G = I.C' * I.A.H
    N = size(G, 2) - 1
    K = eltype(G)

    # forward
    u = fill!(similar(G, size(G, 1), N), zero(K))
	u[:, 1] = α(I.B)
	for n in 2:N
        tmp = u[:, n - 1] .* G[:, n - 1]
		u[:, n] = T(I.B)' * tmp
	end

    # backward
	v = fill!(similar(G, size(G, 1), N), zero(K))
    v[:, N] = G[:, N + 1] .* ω(I.B)
	for n in N - 1: -1: 1
        tmp = G[:, n + 1] .* v[:, n + 1]
		v[:, n] = T(I.B) * tmp
	end

    u, v
end

function Base.sum(I::IntersectedDenseFST, N=size(I.A.H, 2))
    G = I.C' * I.A.H
    u = G[:, 1] .* α(I.B)
    for n in 2:N - 1
        u = (T(I.B)' * u) .* G[:, n]
    end
    dot(u, ω(I.B) .* G[:, end]) + ρ(I)
end

