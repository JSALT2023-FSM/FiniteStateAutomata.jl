# SPDX-License-Identifier: CECILL-2.1

# TODO performance can be improved by preallocating buffers for
# `iterate`.

Base.iterate(A::TransducerOrAcceptor) = nstates(A) == 0 ? nothing : (α(A), (1, α(A)))

function Base.iterate(A::TransducerOrAcceptor, state)
    n, uₙ = state
    n > nstates(A) && throw(ArgumentError("Iterating over a cyclic FST"))
    uₙ₊₁ = transpose(transpose(uₙ) * T(A))
    nnz(uₙ₊₁) > 0 ? (uₙ₊₁, (n+1, uₙ₊₁)) : nothing
end

"""
    W(A)

Return the total weight of `A`, i.e. the ``\\oplus``-sum of all the
path's weight in `A`.
"""
function W(A::TransducerOrAcceptor)
    # ρ(A) ⊕ dot(α(A), MatrixPowerSeries(T(A)), ω(A))
    W = ρ(A)
    for uₙ in A
        W += dot(uₙ, ω(A))
    end
    W
end

function ChainRulesCore.rrule(::typeof(Base.sum), A::AbstractFST{LogSemiring{Tw}}) where Tw
    K = LogSemiring{Tw}

    W = ρ(A)
    U = []
    for uₙ in A
        push!(U, uₙ)
        W += dot(uₙ, ω(A))
    end

    function pullback(Δy)
        V = []
        for vₙ in reverse(A)
            push!(V, vₙ)
        end

        ū = sum(U)
        v̄ = sum(V)
        I_α, _ = findnz(α(A))
        I_ω, _ = findnz(ω(A))

        I_T, J_T, _ = findnz(T(A))
        V_T = zeros(K, length(J_T))
        k_1 = length(U)
        for i in 1:k_1
            for j in 1:(k_1 - i)
                uᵢ, vⱼ = U[i], V[j]
                V_T .+= uᵢ[I_T] .* vⱼ[J_T]
            end
        end

        f = exp ∘ val
        Δx = Δy / f(W)
        ΔA = FST(
            sparsevec(I_α, Δx * f.(v̄[I_α]), nstates(A)),
            sparse(I_T, J_T, Δx * f.(V_T), nstates(A), nstates(A)),
            sparsevec(I_ω, Δx * f.(ū[I_ω]), nstates(A)),
            iszero(ρ(A)) ? zero(Tw) : Δx,
            λ(A)
        )

        (NoTangent(), ΔA)
    end

    W, pullback
end

function ChainRulesCore.rrule(::typeof(Base.sum), A::AbstractFST{TropicalSemiring{Tw}}) where Tw
    W = ρ(A)
    U = []
    for uₙ in A
        push!(U, uₙ)
        W += dot(uₙ, ω(A))
    end

    function pullback(Δy)
        K = BoolSemiring

        V = []
        for vₙ in reverse(A)
            push!(V, vₙ)
        end

        ū = sum(U)
        v̄ = sum(V)
        I_α, _ = findnz(α(A))
        I_ω, _ = findnz(ω(A))

        I_T, J_T, _ = findnz(T(A))
        V_T = zeros(K, length(J_T))
        k_1 = length(U)
        for i in 1:k_1
            for j in 1:(k_1 - i)
                uᵢ, vⱼ = U[i], V[j]
                V_T .+= uᵢ[I_T] .* vⱼ[J_T]
            end
        end

        f = exp ∘ val
        Δx = Δy / f(W)
        ΔA = FST(
            sparsevec(I_α, Δx * f.(v̄[I_α]), nstates(A)),
            sparse(I_T, J_T, Δx * f.(V_T), nstates(A), nstates(A)),
            sparsevec(I_ω, Δx * f.(ū[I_ω]), nstates(A)),
            iszero(ρ(A)) ? zero(Tw) : Δx,
            λ(A)
        )

        (NoTangent(), ΔA)
    end

    W, pullback
end

