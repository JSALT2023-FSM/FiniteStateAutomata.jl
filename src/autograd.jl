# SPDX-License-Identifier: CECILL-2.1
"""
Computes gradient of sum(denseFST ∩ FST) w.r.t DenseFST
"""
function gradient(::typeof(sum), I::IntersectedDenseFST)
    u, v = forward_backward(I)
    N = size(u, 2)
    W = dot(u[:, N] .* (I.C' * I.A.H[:, N]), v[:, N]) + ρ(I)

    ∇ₕW = fill!(similar(I.A.H), zero(eltype(I.A.H)))
    ∇ₕW[:, 1:end-1] = I.C * (u .* v)  # ∂h1, …, ∂hN
    ∇ₕW[:, end] = I.C * ((I.C' * I.A.H[:, end - 1]) .* u[:, end] .* ω(I.B))  # ∂hN+1
    return W, ∇ₕW
end
