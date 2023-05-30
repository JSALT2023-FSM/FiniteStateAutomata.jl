# SPDX-License-Identifier: CECILL-2.1

"""
    W(A)

Return the total weight of `A`, i.e. the ``\\oplus``-sum of all the
path's weight in `A`.
"""
W(A::TransducerOrAcceptor) = ρ(A) ⊕ dot(α(A), MatrixPowerSum(T(A)), ω(A))

