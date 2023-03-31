struct DenseFSA{K, L} <: AbstractAcyclicFSA{K, L}
    H::AbstractMatrix{K}
    Σ::AbstractVector{L}
    ρ::K
    DenseFSA(H, Σ, ρ) = length(Σ) == size(H, 1) ? new{eltype(H), eltype(Σ)}(H, Σ, ρ) : error("Σ not compatible")
end

α(G::DenseFSA) = begin
    L = length(G.Σ)
    N = size(G.H, 2)
    S = L * (N-1)
    sparsevec(1:L, G.H[:, 1], S)
end
T(G::DenseFSA{K}) where K = begin
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
ω(G::DenseFSA) = begin
    L = length(G.Σ)
    N = size(G.H, 2)
    S = L*(N - 1)
    sparsevec(S - L + 1:S, G.H[:, N], S)
end
ρ(G::DenseFSA) = G.ρ
λ(G::DenseFSA) = repeat(G.Σ, size(G.H, 2) - 1)

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
    # Hn = C' * A.H
    # Hn[:, 1] .*= α(B)
    # Hn[:, end] .*= ω(B)
    # DenseFSA(Hn, λ(B), ρ(A) * ρ(B))
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
