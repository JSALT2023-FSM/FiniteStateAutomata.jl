struct DenseFST{K, L} <: AbstractAcyclicFST{K, L}
    H::AbstractMatrix{K}
    Σ::AbstractVector{L}
    ρ::K
    DenseFST(H, Σ, ρ) = length(Σ) == size(H, 1) ? new{eltype(H), eltype(Σ)}(H, Σ, ρ) : error("Σ not compatible")
end

DenseFST(H::AbstractMatrix{K}, Σ::AbstractVector) where K = DenseFST(H, Σ, zero(K))

α(G::DenseFST) = begin
    L = length(G.Σ)
    N = size(G.H, 2)
    S = L * (N-1)
    sparsevec(1:L, G.H[:, 1], S)
end
T(G::DenseFST{K}) where K = begin
    rows = []
    cols = []
    vals = K[]
    L = length(G.Σ)
    N = size(G.H, 2)
    S = L * (N - 1)
    for n in 2:N - 1
        v = repeat(G.H[:, n], L) # flatten
        r = [(n-2)*L + i for i in 1:L for _ in 1:L]
        c = [(n-1)*L + i for _ in 1:L for i in 1:L]
        rows = vcat(rows, r)
        cols = vcat(cols, c)
        vals = vcat(vals, v)
    end
    sparse(rows, cols, vals, S, S)
end
ω(G::DenseFST) = begin
    L = length(G.Σ)
    N = size(G.H, 2)
    S = L*(N - 1)
    sparsevec(S - L + 1:S, G.H[:, N], S)
end
ρ(G::DenseFST) = G.ρ
λ(G::DenseFST) = repeat(G.Σ, size(G.H, 2) - 1)

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
    # Hn = C' * A.H
    # Hn[:, 1] .*= α(B)
    # Hn[:, end] .*= ω(B)
    # DenseFST(Hn, λ(B), ρ(A) * ρ(B))
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

function Base.sum(I::IntersectedDenseFST, n = size(I.A.H, 2) - 1)
    CH = I.C' * I.A.H
    v = CH[:, 1] .* α(I.B)
    for n in 2:n
        v = (T(I.B)' * v) .* CH[:, n]
    end
    dot(v, ω(I.B) .* CH[:, end]) + ρ(I)
end

