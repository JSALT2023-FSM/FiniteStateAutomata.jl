# SPDX-License-Identifier: CECILL-2.1

Base.iterate(A::AbstractFST) = nstates(A) == 0 ? nothing : (α(A), α(A))

function Base.iterate(A::AbstractFST, uₙ)
    uₙ₊₁ = T(A)' * uₙ
    nnz(uₙ₊₁) > 0 ? (uₙ₊₁, uₙ₊₁) : nothing
end

