# SPDX-License-Identifier: CECILL-2.1
#
using CUDA
using FiniteStateAutomata
using LinearAlgebra
using Random
using BenchmarkTools
using Random

include("abstracttypes.jl")

import FiniteStateAutomata.nzrange
import FiniteStateAutomata.colvals
import FiniteStateAutomata.getrowptr
import FiniteStateAutomata.nonzeros
import FiniteStateAutomata.nnz

struct CuSparseMatrixCSR{S} <: AbstractSparseMatrixCSR{S}
    m::Int
    n::Int
    rowptr::CuArray{Int, 1}
    colval::CuArray{Int, 1}
    nzval::CuArray{S, 1}
end

function to_cpu(x::CuSparseMatrixCSR{S}) where S
    return SparseMatrixCSR(x.m, x.n
                             Array(x.rowptr), 
                             Array(x.colval), 
                             Array(x.nzval))
end

function to_gpu(x::SparseMatrixCSR{S}) where S
    return CuSparseMatrixCSR(x.m, x.n
                             CUDA.CuArray(x.rowptr), 
                             CUDA.CuArray(x.colval), 
                             CUDA.CuArray(x.nzval))
end

nnz(v::CuSparseMatrixCSR) = length(v.nzval)


function Base.size(x::CuSparseMatrixCSR{S}) where S
    return (x.m, x.n)
end

CuSparseMatrixCSR(A::SparseMatrixCSR) = CuSparseMatrixCSR(
                                                   A.m,
                                                   A.n,
                                                   CuArray(A.rowptr),
                                                   CuArray(A.colval),
                                                   CuArray(A.nzval)
                                                  )
function nzrange(X::CuSparseMatrixCSR, r) 
    q = Array(X.rowptr)
    return  q[r]:(q[r+1]-1)
end

getrowptr(X::CuSparseMatrixCSR) = X.rowptr
colvals(X::CuSparseMatrixCSR) = X.colval
nonzeros(X::CuSparseMatrixCSR) = X.nzval


CUDA.@allowscalar function Base.getindex(X::CuSparseMatrixCSR{S}, i::Int, j::Int) where S
    ii = nzrange(X, i)
    isempty(ii) && return zero(S)
    nnzy = colvals(X)[ii]
    fj = findfirst(a -> a==j, nnzy)
    (isnothing(fj)) ? zero(S) : nonzeros(X)[ii[1]-1+fj]
end

struct CuSparseVectorX{S} <: AbstractSparseVector{S}
    n::Int
    nzind::CuArray{Int, 1}
    nzval::CuArray{S, 1}
end


function Base.size(x::CuSparseVectorX{S}) where S
    return (x.n, )
end


CuSparseVectorX(A::SparseVector) = CuSparseVectorX(
                                               n(A.n),
                                               CuArray(A.nzind),
                                               CuArray(A.nzval)
                                              )

# Create an uninitialized sparse vector with internal buffers to
# store `nzc` non-zero values.
CuSparseVectorX(S, n, nzc) =
    CuSparseVectorX(n, CuArray{Int}(undef, nzc), CuArray{S}(undef, nzc))

function CuSparseVectorX(xt::LinearAlgebra.Transpose{S, <:SparseVector}) where S
    x = parent(xt)
    y = CuSparseVectorX(x.n, CuArray(x.nzind), CuArray(x.nzval))
    return LinearAlgebra.Transpose(y)
end

function Base.getindex(x::CuSparseVectorX{S}, i::Int) where S
    nzi = findfirst(a -> a == i, x.nzind)
    isnothing(nzi) ? zero(S) : x.nzval[nzi]
end

         
nnz(v::CuSparseVectorX) = length(v.nzval)

function _k_cu_spmv!(m, n,
        xnzind, xnzval,
        arowptr, acolval, anzval,
        yout_dense::CuDeviceVector{S, 1}
    ) where S
   
    r = blockIdx().x # each block takes care of a single row
    t = threadIdx().x 
    nof_threads = blockDim().x
    idx = t + (r - 1)*nof_threads

    if r <= n
        yout_dense[r] = 0
    end

    if idx > length(xnzind)
        return
    end

    rowx = xnzind[idx]
    col_idx = arowptr[rowx]
    col_idxp1 = arowptr[rowx + 1] - 1
    for i in col_idx:col_idxp1
        CUDA.@atomic yout_dense[acolval[i]] = yout_dense[acolval[i]]  ⊕ (xnzval[idx] ⊗ anzval[i])
    end

    #offset = warpsize()/2
    #while offset > 0
    #     accum = accum ⊕ shfl_down_sync(CUDA.FULL_MASK, accum,  offset)
    #     offset = offset/2 
    #end
    #yout_dense[r] = accum
    return
end 

# more optimized function (the warp reads continuous memoryin col_idx)
@inbounds function _k_cu_spmv2!(m, n,
        xnzind, xnzval,
        arowptr, acolval, anzval,
        yout_dense::CuDeviceVector{S, 1}
    ) where S
   
    r = blockIdx().x # each block takes care of a single row
    t = threadIdx().x 
    nof_threads = blockDim().x

    tid = t + (r - 1) * nof_threads

    if tid <= n
        yout_dense[tid] = 0
    end

    rowx = xnzind[r]
    col_idx = arowptr[rowx] + t - 1
    col_idxp1 = arowptr[rowx + 1] - 1
    for i in col_idx:nof_threads:col_idxp1
        CUDA.@atomic yout_dense[acolval[i]] = yout_dense[acolval[i]]  ⊕ (xnzval[r] ⊗ anzval[i])
    end

    return
end 

function mult_spvspm!(y_out::CuArray{S}, xᵀ::Transpose{S,<:CuSparseVectorX}, A::CuSparseMatrixCSR{S}) where S
    size(xᵀ, 2) != size(A, 1) && throw(DimensionMismatch())
    x = parent(xᵀ)
    kernel = @cuda launch=false _k_cu_spmv!(A.m, A.n, x.nzind, x.nzval,
                                          A.rowptr, A.colval, A.nzval, y_out)
   
    config = launch_configuration(kernel.fun)

    @assert length(x.nzval) < (2^31)
    blocks = length(x.nzval) #internal limit 2^31 - 1
    threads = config.threads
    if (A.n < threads)
        threads = A.n
    end

    kernel(A.m, A.n, x.nzind, x.nzval,
          A.rowptr, A.colval, A.nzval, y_out; 
          threads=threads, blocks=blocks)
    nzind = findall(!iszero, y_out)
    nzval = y_out[nzind]
    
    out = CuSparseVectorX(A.n, nzind, nzval)
    synchronize()    

    return out
end

function Base.:*(xᵀ::Transpose{S,<:CuSparseVectorX}, A::CuSparseMatrixCSR{S}) where S
    y_out = CuArray(Array{S}(undef, A.n))
    out = mult_spvspm!(y_out, xᵀ, A)
    return out
end



CUDA.allowscalar(true)
rng = MersenneTwister(1234)

A = Float32.(rand(rng, [-1, 0, 0, 0], (10000, 10000)))
x = Float32.(rand(rng, [-1, 0, 0, 0, 0,], (1, 10000))) 

#A = Float32[-1.0 -2.0; 3.0 -4.0; 5.0 -6.0]
#x = transpose(Float32[-1, 2, -3])
#@show A
#@show x

xs = sparsevec(x)
As = sparse(A)

println("CPU")
#@show (xs, parent(xs).nzind, parent(xs).nzval)
#@show As, As.rowptr, As.colval, As.nzval
@btime xs * As

cuxs = CuSparseVectorX(xs)
cuas = CuSparseMatrixCSR(As)

#@show typeof(cuxs), cuxsI
#@show typeof(cuas), cuas

cuys = cuxs * cuas
println("CUDA")
@btime cuxs * cuas
#@show typeof(cuys), cuys

println("CUDA with buffer")
y_out = CuArray(Array{Float32}(undef, As.n))
mult_spvspm!(y_out, cuxs, cuas)
@btime mult_spvspm!(y_out, cuxs, cuas)
#@show typeof(cuys), cuys
#o = reshape(Array(Float32.(cuys)), 1, :)
#p = Float32.(Array(x*A))

#@show o
#@show p
#@assert o == p
