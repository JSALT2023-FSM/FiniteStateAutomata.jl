# SPDX-License-Identifier: CECILL-2.1

# Generic intersection
struct IntersectedFSA{K, L, T_A<:AbstractFSA{K,L}, T_B<:AbstractFSA{K,L}} <: AbstractFSA{K, L}
    A::T_A
    B::T_B
    M::AbstractMatrix{K}
    λ::AbstractVector{L}
end

_computeM(::Type{K}, λ₁::AbstractVector{L}, λ₂::AbstractVector{L}) where {K,L} = begin
    rows = []
    labels = L[]
    q = 0
    for (i, s1) in enumerate(λ₁)
        for (j, s2) in enumerate(λ₂)
            if s1 == s2
                q += 1
                push!(rows, (i-1) * length(λ₂) + j)
                push!(labels, s1)
            end
        end
    end
    sparse(rows, 1:q, one(K), length(λ₁) * length(λ₂), q), labels
end

IntersectedFSA(A::AbstractFSA{K}, B::AbstractFSA{K}) where K = begin 
    M, labels = _computeM(K, λ(A), λ(B))
    IntersectedFSA(A, B, M, labels)
end

Base.intersect(A::AbstractFSA, B::AbstractFSA) = IntersectedFSA(A, B)

α(I::IntersectedFSA) = I.M' * kron(α(I.A), α(I.B))
T(I::IntersectedFSA) = I.M' * kron(T(I.A), T(I.B)) * I.M
ω(I::IntersectedFSA) = I.M' * kron(ω(I.A), ω(I.B))
ρ(I::IntersectedFSA) = ρ(I.A) * ρ(I.B)

# Intersection with DenseFSA
struct IntersectedDenseFSA{K, L, T<:AbstractFSA{K, L}} <: AbstractAcyclicFSA{K, L}
    A::DenseFSA{K, L}
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

IntersectedDenseFSA(A::DenseFSA{K}, B::AbstractFSA{K}) where K = begin
    C = _computeC(K, A.Σ, λ(B)) # TODO switch cols and rows
    IntersectedDenseFSA(A, B, C)
end

Base.intersect(A::DenseFSA, B::AbstractFSA) = IntersectedDenseFSA(A, B) 
Base.intersect(A::AbstractFSA, B::DenseFSA) = IntersectedDenseFSA(B, A) # we assume ∩ to be commutative w.r.t any K

α(I::IntersectedDenseFSA) = begin
    h1 = I.A.H[:, 1]
    Q = nstates(I.B) * (size(I.A.H, 2) - 1)
    vals = (I.C' * h1) .* α(I.B)
    sparsevec(findnz(vals)..., Q)
end

ω(I::IntersectedDenseFSA) = begin
    h = I.A.H[:, end]
    Q = nstates(I.B) * (size(I.A.H, 2) - 1)
    vals = (I.C' * h) .* ω(I.B)
    I, V = findnz(vals)
    sparsevec(I .+ (Q - length(vals)), V, Q)
end

T(I::IntersectedDenseFSA) = begin
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

λ(I::IntersectedDenseFSA) = repeat(λ(I.B), size(I.A.H, 2) - 1)
ρ(I::IntersectedDenseFSA) = ρ(I.B) * ρ(I.A)

function forward_backward(I::IntersectedDenseFSA)
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

function Base.sum(I::IntersectedDenseFSA, N=size(I.A.H, 2))
    G = I.C' * I.A.H
    u = G[:, 1] .* α(I.B)
    for n in 2:N - 1
        u = (T(I.B)' * u) .* G[:, n]
    end
    dot(u, ω(I.B) .* G[:, end]) + ρ(I)
end
