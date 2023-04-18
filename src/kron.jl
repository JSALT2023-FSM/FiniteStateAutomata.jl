# SPDX-License-Identifier: CECILL-2.1

struct KronFST{K, L, T_A<:AbstractFST{K,L}, T_B<:AbstractFST{K,L}} <: AbstractFST{K, L}
    f
    A::T_A
    B::T_B
end

Base.kron(f, A::AbstractFST, B::AbstractFST) = KronFST(f, A, B)
Base.kron(A::AbstractFST, B::AbstractFST) = kron((x, y) -> "($x, $y)", A, B)

function ChainRulesCore.rrule(::typeof(Base.kron), f, A::AbstractFST{K}, B::AbstractFST{K}) where K
    kC = KronFST(f, A, B)

    function pullback(ΔY)
        N = nstates(A)
        M = nstates(B)

        I_A, V_A_old = findnz(α(A))
        I_B, V_B_old = findnz(α(B))
        V_A = fill!(similar(V_A_old), zero(K))
        V_B = fill!(similar(V_B_old), zero(K))
        for n in 1:length(I_A)
            for m in 1:length(I_B)
                i = (I_A[n]-1) * M + I_B[m]
                V_A[n] += α(ΔY)[i] * V_B_old[m]
                V_B[m] += α(ΔY)[i] * V_A_old[n]
            end
        end
        α_A = sparsevec(I_A, V_A, N)
        α_B = sparsevec(I_B, V_B, M)

        I_A, J_A, V_A_old = findnz(T(A))
        I_B, J_B, V_B_old = findnz(T(B))
        V_A = fill!(similar(V_A_old), zero(K))
        V_B = fill!(similar(V_B_old), zero(K))
        for n in 1:length(I_A)
            for m in 1:length(I_B)
                i = (I_A[n]-1) * M + I_B[m]
                j = (J_A[n]-1) * M + J_B[m]
                V_A[n] += T(ΔY)[i,j] * V_B_old[m]
                V_B[m] += T(ΔY)[i,j] * V_A_old[n]
            end
        end
        T_A = sparse(I_A, J_A, V_A, N, N)
        T_B = sparse(I_B, J_B, V_B, M, M)

        I_A, V_A_old = findnz(ω(A))
        I_B, V_B_old = findnz(ω(B))
        V_A = fill!(similar(V_A_old), zero(K))
        V_B = fill!(similar(V_B_old), zero(K))
        for n in 1:length(I_A)
            for m in 1:length(I_B)
                i = (I_A[n]-1) * M + I_B[m]
                V_A[n] += ω(ΔY)[i] * V_B_old[m]
                V_B[m] += ω(ΔY)[i] * V_A_old[n]
            end
        end
        ω_A = sparsevec(I_A, V_A, N)
        ω_B = sparsevec(I_B, V_B, M)

        ΔA = FST(α_A, T_A, ω_A, ρ(B) * ρ(ΔY), λ(A))
        ΔB = FST(α_B, T_B, ω_B, ρ(A) * ρ(ΔY), λ(B))

        (NoTangent(), NoTangent(), ΔA, ΔB)
    end

    kC, pullback
end

α(kA::KronFST) = kron(α(kA.A), α(kA.B))
T(kA::KronFST) = kron(T(kA.A), T(kA.B))
ω(kA::KronFST) = kron(ω(kA.A), ω(kA.B))
ρ(kA::KronFST) = ρ(kA.A) * ρ(kA.B)
λ(kA::KronFST) = vec([kA.f(l_a, l_b) for l_b in λ(kA.B), l_a in λ(kA.A)])

