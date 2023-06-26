# SPDX-License-Identifier: CECILL-2.1

struct SparseVector{S} <: AbstractSparseVector{S}
    n::Int
    nzind::Vector{Int}
    nzval::Vector{S}
end

# Create an uninitialized sparse vector with internal buffers to
# store `nzc` non-zero values.
SparseVector(S::Type{<:Semiring}, n, nzc) =
    SparseVector(n, Vector{Int}(undef, nzc), Vector{S}(undef, nzc))

function gather!(x::SparseVector, acc::SparseAccumulator)
    length(x.nzind) != nnz(acc) && throw(DimensionMistmatch())

    for n in 1:nnz(acc)
        i = acc.ind[n]
        x.nzind[n] = i
        x.nzval[n] = acc.val[i]
    end

    x
end

#= AbstractArray API =#

Base.size(x::SparseVector) = (x.n,)

function Base.getindex(x::SparseVector{S}, i::Int) where S
    nzi = findfirst(a -> a == i, x.nzind)
    isnothing(nzi) ? zero(S) : x.nzval[nzi]
end

#= Sparse API =#

rowvals(x::SparseVector) = x.nzind
nonzeros(x::SparseVector) = x.nzval

function sparsevec(I::AbstractVector, V::AbstractVector{S}, n = maximum(I)) where S
    acc = SparseAccumulator(eltype(V), n)
    for i in 1:length(I) scatter!(acc, V[i], I[i]) end
    spv = SparseVector(S, n, nnz(acc))
    gather!(spv, acc)
end

sparsevec(I::AbstractVector, val, n = maximum(I)) = sparsevec(I, repeat([val], n), n)

nnz(x::SparseVector) = length(x.nzind)
nzrange(x::SparseVector) = 1:nnz(x)
findnz(x::SparseVector) = (x.nzind, x.nzval)

