# SPDX-License-Identifier: CECILL-2.1

"""
    Base.sum(A; n = nstates(A))

Accumulate the weights of all the paths. For cyclic FSTs sum after `n`
steps.
"""
function Base.sum(A::AbstractFST; n = nstates(A))
    W = ρ(A)
    for (i, uₙ) in enumerate(A)
        i <= n || break
        W += dot(uₙ, ω(A))
    end
    W
end

function ChainRulesCore.rrule(::typeof(Base.sum), A::AbstractFST{K}) where K
    W = ρ(A)
    U = []
    for uₙ in A
        push!(U, uₙ)
        W += dot(uₙ, ω(A))
    end

    function pullback(ΔY)
        V = []
        for vₙ in reverse(A)
            push!(V, vₙ)
        end

        ū = sum(U)
        I_α, _ = findnz(α(A))

        I_T, J_T, _ = findnz(T(A))
        V_T = zeros(K, length(J_T))
        k_1 = length(U)
        for i in 1:k_1
            for j in 1:(k_1 - i)
                uᵢ, vⱼ = U[i], V[k_1 - j]
                V_T .+= uᵢ[I_T] .* vⱼ[J_T]
            end
        end

        v̄ = sum(V)
        I_ω, _ = findnz(ω(A))

        ΔA = FST(
            sparsevec(I_α, ū[I_α], nstates(A)),
            sparse(I_T, J_T, V_T, nstates(A), nstates(A)),
            sparsevec(I_ω, v̄[I_ω], nstates(A)),
            iszero(ρ(A)) ? zero(K) : one(K),
            λ(A)
        )

        (NoTangent(), ΔA)
    end

    W, pullback
end

