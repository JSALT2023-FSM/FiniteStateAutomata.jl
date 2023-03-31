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

