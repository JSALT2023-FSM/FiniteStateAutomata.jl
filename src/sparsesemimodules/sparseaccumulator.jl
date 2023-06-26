# SPDX-License-Identifier: CECILL-2.1

mutable struct SparseAccumulator{TV<:AbstractVector,
								 TM<:AbstractArray{Bool},
								 TI<:AbstractArray{<:Integer},
								 Ti}
	val::TV
	mask::TM
	ind::TI
	nnz::Ti
end

# Create a sparse accumulator of type `T` with internal buffer of size
# `n`. `m` is the maximum number of nonzero values.
SparseAccumulator(T, n::Integer, m = n) = SparseAccumulator(
    Vector{T}(undef, n),
    zeros(Bool, n),
    Vector{Int}(undef, m),
    Threads.Atomic{Int}(0)
)

function scatter!(acc::SparseAccumulator, v, i::Integer)
	if ! acc.mask[i]
		acc.mask[i] = true
		acc.val[i] = v
		nzi = Threads.atomic_add!(acc.nnz, 1)
		acc.ind[nzi+1] = i
	else
		acc.val[i] = acc.val[i] ⊕ v
	end
end

function reset!(acc::SparseAccumulator)
    acc.mask[acc.ind[1:nnz(acc)]] .= false
    acc.nnz = Threads.Atomic{Int}(0)
end

nnz(acc::SparseAccumulator) = acc.nnz[]

