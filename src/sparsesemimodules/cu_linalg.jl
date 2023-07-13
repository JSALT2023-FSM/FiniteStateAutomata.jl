# SPDX-License-Identifier: CECILL-2.1
#

struct CuSparseMatrixCSR{S} <: AbstractSparseMatrixCSR{S}
    m::Int
    n::Int
    rowptr::CuArray{Int, 1}
    colval::CuArray{Int, 1}
    nzval::CuArray{S, 1}
end

function to_cpu(x::CuSparseMatrixCSR{S}) where S
    return SparseMatrixCSR(x.m, x.n,
                             Array(x.rowptr),
                             Array(x.colval),
                             Array(x.nzval))
end

function to_gpu(x::SparseMatrixCSR{S}) where S
    return CuSparseMatrixCSR(x.m, x.n,
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

function to_cpu(x::CuSparseVectorX{S}) where S
    return SparseVector(x.n,
                     Array(x.nzind),
                     Array(x.nzval))
end

function to_cpu(xᵀ::Transpose{S, <:CuSparseVectorX{S}}) where S
    x = parent(xᵀ)
    return transpose(SparseVector(x.n,
                             Array(x.nzind),
                             Array(x.nzval)))
end

function to_gpu(x::SparseVector{S}) where S
    return CuSparseVectorX(x.n,
                         CUDA.CuArray(x.nzind),
                         CUDA.CuArray(x.nzval))
end

function to_gpu(xᵀ::Transpose{S, <:SparseVector{S}}) where S
    x = parent(xᵀ)
    return transpose(CuSparseVectorX(x.n,
                         CUDA.CuArray(x.nzind),
                         CUDA.CuArray(x.nzval)))
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

    if idx <= length(yout_dense)
        yout_dense[idx] = 0
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
function _k_cu_spmv2!(m, n,
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

    #Completely zero matrix
    if length(x.nzval) == 0
        nzval = CUDA.cu(S[])
        nzind = CUDA.cu(Int[])
        out = CuSparseVectorX(A.n, nzind, nzval)
        return out
    end

    kernel = @cuda launch=false _k_cu_spmv2!(A.m, A.n, x.nzind, x.nzval,
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

    return transpose(out)
end

@inbounds function _k_cu_spmspm!(a_m, a_n, row, arowptr, acolval, anzval,
                       b_m, b_n, browptr, bcolval, bnzval,
                       y_buffer)

    r = blockIdx().x # each block takes care of a single row
    t = threadIdx().x
    nof_threads = blockDim().x
    tid = t + (r - 1)*nof_threads

    if tid <= length(y_buffer)
        y_buffer[tid] = 0
    end

    rangeptr = arowptr[row]:(arowptr[row+1]-1)
    xnzind = view(acolval, rangeptr)
    xnzval = view(anzval, rangeptr)

    if  r > length(xnzind)
        return
    end

    rowx = xnzind[r]
    col_idx = browptr[rowx] + t - 1
    col_idxp1 = browptr[rowx + 1] - 1
    for i in col_idx:nof_threads:col_idxp1
        q =  bcolval[i]
        CUDA.@atomic y_buffer[q] = y_buffer[q]  ⊕ (xnzval[r] ⊗ bnzval[i])
    end

    return

end

@inbounds function _k_cu_spmspm3!(a_m, a_n, row, arowptr, acolval, anzval,
                       b_m, b_n, browptr, bcolval, bnzval,
                       y_buffer)

    r = blockIdx().x # each block takes care of a single row
    t = threadIdx().x
    nof_threads = blockDim().x
    tid = t + (r - 1)*nof_threads

    if tid <= length(y_buffer)
        y_buffer[tid] = 0
    end

    rangeptr = arowptr[row]:(arowptr[row+1]-1)
    xnzind = view(acolval, rangeptr)
    xnzval = view(anzval, rangeptr)

    if  r > length(xnzind)
        return
    end

    rowx = xnzind[r]
    col_idx = browptr[rowx] + t - 1
    col_idxp1 = browptr[rowx + 1] - 1
    for i in col_idx:nof_threads:col_idxp1
        q =  bcolval[i]
        CUDA.@atomic y_buffer[q] = y_buffer[q]  ⊕ (xnzval[r] ⊗ bnzval[i])
    end

    return

end

function _k_assign!(nzval::CUDA.CuDeviceVector{Float32, 1},
        colval::CUDA.CuDeviceVector{Int64, 1},
        rowptr::CUDA.CuDeviceVector{Int64, 1},
        row::Int,
        row_offset::Int,
        nzind::CUDA.CuDeviceVector{Int, 1},
        y_buffer::CUDA.CuDeviceVector{Float32, 1})
    r = blockIdx().x # each block takes care of a single row
    t = threadIdx().x
    nof_threads = blockDim().x
    k = length(nzind)

    tid = t + (r - 1) * nof_threads
    if tid <= k
        v = nzind[tid]

        nzval[row_offset + tid] = y_buffer[v]
        colval[row_offset + tid] = v
    end
    if tid == 1
        rowptr[row]  = row_offset  + 1
        rowptr[row+1]  = row_offset + k +1
    end
    return
end

function mult_spmspm!(y_buffer::CuArray{S}, A::CuSparseMatrixCSR{S}, B::CuSparseMatrixCSR{S}) where S
    size(A, 2) != size(B, 1) && throw(DimensionMismatch())

    #Completely zero matrix
    nnza = nnz(A)
    nnzb = nnz(B)
    max_elems = min(nnza * nnzb, A.m * B.n)
    #@show nnza, nnzb
    rowptr = CUDA.zeros(Int, A.m+1)
    colval = CUDA.zeros(Int, max_elems)
    nzval = CUDA.zeros(S, max_elems)

    #rowptr = CuArray(Array{Int}(undef, A.m))
    #colval = CuArray(Array{Int}(undef, max_elems))
    #nzval = CuArray(Array{S}(undef, max_elems))

    if nnza == 0 || nnzb == 0
        out = CuSparseMatrixCSR(A.m, B.n, rowptr, colval, nzval)
        return out
    end

    row_offset = 0
    for row in eachindex(A.rowptr[1:end-1])
        kernel = @cuda  launch=false  _k_cu_spmspm!(
                                             A.m, A.n, row,
                                             A.rowptr, A.colval, A.nzval,
                                             B.m, B.n,
                                             B.rowptr, B.colval, B.nzval,
                                             y_buffer)
        config = launch_configuration(kernel.fun)

        blocks = B.m #internal limit 2^31 - 1
        threads = config.threads
        if (B.n < threads)
            threads = B.n
        end
        kernel(
             A.m, A.n, row,
             A.rowptr, A.colval, A.nzval,
             B.m, B.n,
             B.rowptr, B.colval, B.nzval,
             y_buffer; blocks=blocks, threads=threads)

        out_nzind = findall(!iszero, y_buffer)
        #if length(out_nzind) == 0
        #    continue
        #end
        threads = 512
        blocks = cld(length(out_nzind)+1, threads)
        @cuda threads=threads blocks=blocks _k_assign!(nzval, colval, rowptr,
                                             row, row_offset,
                                             out_nzind, y_buffer)
        row_offset += length(out_nzind)
    end
    return CuSparseMatrixCSR(A.m, B.n, rowptr[1:(A.m+1)], colval[1:row_offset], nzval[1:row_offset])
end


function Base.:*(xᵀ::Transpose{S,<:CuSparseVectorX}, A::CuSparseMatrixCSR{S}) where S
    y_out = CuArray(Array{S}(undef, A.n))
    out = mult_spvspm!(y_out, xᵀ, A)
    return out
end

function Base.:*(A::CuSparseMatrixCSR{S}, B::CuSparseMatrixCSR{S}) where S
    y_out = CuArray(Array{S}(undef, B.n))
    out = mult_spmspm!(y_out, A, B)
    return out
end

