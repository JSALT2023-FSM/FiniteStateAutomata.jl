struct DenseFSA{K, L} <: AbstractAcyclicFSA{K, L}
    H::AbstractMatrix{K}
    Σ::AbstractVector{L}
    ρ::K
    DenseFSA(H, Σ, ρ) = length(Σ) == size(H, 1) ? new{eltype(H), eltype(Σ)}(H, Σ, ρ) : error("Σ not compatible")
end

DenseFSA(H::AbstractMatrix{K}, Σ::AbstractVector) where K = DenseFSA(H, Σ, zero(K))

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

