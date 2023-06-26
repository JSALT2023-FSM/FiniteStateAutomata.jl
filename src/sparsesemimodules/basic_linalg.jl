# SPDX-License-Identifier: CECILL-2.1

#=

Low-level functions for basic linear algebra operations for CPU and
GPU.

=#


#= Scalar OP =#

function _apply_scalar_op!(y::AbstractArray, op::Function, x::AbstractArray, a)
    length(y) == length(x) || throw(DimensionMismatch())
    for i in eachindex(x)
        y[i] = op( x[i], a )
    end
    y
end

function _apply_scalar_op!(y::CuArray, op::Function, x::CuArray, a)
    length(y) == length(x) || throw(DimensionMismatch())

    kernel = @cuda launch=false _k_asop!(y, op, x, a)

    threads, blocks = get_config(kernel, length(x))

    CUDA.@sync kernel(y, op, x, a; threads=threads, blocks=blocks)
end

#TODO: Test

function _k_asop!(y, op, x, a)
    idx = threadIdx().x + (blockIdx().x - 1)*blockDim().x

    if idx <= length(x)
		@inbounds y[idx] = op(x[idx], a)
	end

	nothing
end

#= Dot product =#

function _dot(op::Function, x_nzind, x_nzval, y_nzind, y_nzval)
    acc = zero(eltype(x_nzval))
    for i in eachindex(x_nzind)
        xᵢ = x_nzind[i]
        j = searchsorted( y_nzind, xᵢ)
        if(!isempty(j))
            xᵥ = x_nzval[i]
            yᵥ = y_nzval[j][1]
            acc = acc ⊕ op(xᵥ, yᵥ)
        end
    end
    acc
end

function _dot(op::Function, x_nzind::CuArray, x_nzval::CuArray, y_nzind::CuArray, y_nzval::CuArray)
    acc = CuArray(zeros(eltype(x_nzval), 1))

    kernel = @cuda launch=false _k_dot!(acc, op, x_nzind, x_nzval, y_nzind, y_nzval)

    threads, blocks = get_config(kernel, length(x_nzind))

    CUDA.@sync kernel(acc, op, x_nzind, x_nzval, y_nzind, y_nzval; threads=threads, blocks=blocks)

    CUDA.@allowscalar return acc[1]
end

function _k_dot!(acc, op, x_nzind, x_nzval, y_nzind, y_nzval)
    idx = get_idx()

    sum = zero(eltype(acc))

    if idx <= length(x_nzval)
        lane = ((idx - 1) % warpsize()) + 1

        xᵢ = x_nzind[idx]
		j = findfirst(a -> a==xᵢ, y_nzind)

        if (!isnothing(j))
			xᵥ = x_nzval[idx]
            yᵥ = y_nzval[j]

            sum = sum ⊕ op(xᵥ, yᵥ)
        end

        all_reduce!(acc, ⊕, sum, lane)

    end

    nothing
end

#= Applying OP =#

function _apply_op!(z::AbstractArray, op::Function, x::AbstractArray, y::AbstractArray)
    (length(y) == length(x) && length(y) == length(z)) || throw(DimensionMismatch())
    for i in eachindex(z)
        z[i] = op( x[i], y[i] )
    end
    z
end

# TODO: fixme
function _apply_op!(z::AbstractArray, op::Function, x::SparseSMVector, y::SparseSMVector)
    (length(y) == length(x) && length(y) == length(z)) || throw(DimensionMismatch())

    K = eltype(z)
    z = fill!(z, zero(K))

    for i in 1:nnz(x)
        xi = x.nzind[i]
        z[xi] = z[xi] ⊕ x.nzval[i]
    end

    for i in 1:nnz(y)
        yi = y.nzind[i]
        z[yi] = z[yi] ⊕ y.nzval[i]
    end

    z
end

#= Kronecker product =#

function _kron!(z::AbstractVector, x::AbstractVector, y::AbstractVector)
    for i in eachindex(x)
        for j in eachindex(y)
            z[(i-1)*length(y) + j] = x[i] ⊗ y[j]
        end
    end
    z
end

function _kron!(z::CuArray, x::CuArray, y::CuArray)
    kernel = @cuda launch=false _k_kron!(z, x, y)

    threads, blocks = get_config(kernel, length(x)*length(y))

    CUDA.@sync kernel(z, x, y; threads=threads, blocks=blocks)

    return z
end

#TODO: Test

function _k_kron!(z, x, y)
    idx = threadIdx().x + (blockIdx().x - 1)*blockDim().x

    if idx <= length(x) * length(y)
        yᵢ = idx%length(y)
        xᵢ = (idx-yᵢ) ÷ length(y) #source : a code of mine that seemed to work in another project
        #z[(i-1)*length(y) + j] = x[i] ⊗ y[j]
        #(( (idx-idx%l(y)) ÷ l(y) ) -1 ) * l(y) + idx%l(y)
        #oof, isn't it just z[idx] ?
        z[idx] = x[xᵢ] ⊗ y[yᵢ]
    end

end

#= General sequential sparse product =#

# y = A * x
function _spmdv_csc!(y::AbstractVector, op::Function, colptr_A, rowval_A, nzval_A, x)
    y = fill!(y, zero(eltype(y)))

    # Might be a way to inverse the loops
    for i in 1:(length(colptr_A)-1)
        for j in colptr_A[i]:(colptr_A[i+1]-1)
            y[rowval_A[j]] = y[rowval_A[j]] ⊕ op(nzval_A[j], x[i])
        end
    end
    return y
end

function _spmdv_csr!(y::AbstractVector, op::Function, rowptr_A, colval_A, nzval_A, x)
    y = fill!(y, zero(eltype(y)))
    for i in 1:(length(rowptr_A)-1)
        for j in rowptr_A[i]:(rowptr_A[i+1]-1)
            y[i] = y[i] ⊕ op(nzval_A[j], x[colval_A[j]])
        end
    end
    return y
end

function _spmspv_csc!(y::AbstractVector, op::Function, colptr_A, rowval_A, nzval_A, x_nzind, x_nzval)
    fill!(y, zero(eltype(y)))
    for ri in eachindex(x_nzind)
        i = x_nzind[ri]
        vᵢ = x_nzval[ri]
        for cptr in colptr_A[i]:(colptr_A[i+1]-1)
            row = rowval_A[cptr]
            val = nzval_A[cptr]
            y[row] = y[row] ⊕ op(val, vᵢ)
        end
    end
    return y
end

function _spmspv_csr!(y::AbstractVector, op::Function, rowptr_A, colval_A, nzval_A, x_nzind, x_nzval)
    fill!(y, zero(eltype(y)))
    for i in eachindex(rowptr_A[1:end-1])
        rangeptr = rowptr_A[i]:(rowptr_A[i+1]-1)
        y[i] = y[i] ⊕ _dot(op, colval_A[rangeptr], nzval_A[rangeptr], x_nzind, x_nzval)
    end
    return y
end

function _spmspm_csc!(C::AbstractMatrix, op::Function, colptr_A, rowval_A, nzval_A, colptr_B, rowval_B, nzval_B)
    for colB in 1:length(colptr_B)-1
        ptrrange = colptr_B[colB]:(colptr_B[colB+1]-1)
        _spmspv_csc!(view(C, :, colB), op, colptr_A, rowval_A, nzval_A, rowval_B[ptrrange], nzval_B[ptrrange])
    end

    return C
end

function _spmspm_csr!(C::AbstractMatrix, op::Function, rowptr_A, colval_A, nzval_A, rowptr_B, colval_B, nzval_B)
    fill!(C, zero(eltype(C)))

    for rowB in eachindex(rowptr_B[1:end-1])
        for icolB in rowptr_B[rowB]:(rowptr_B[rowB+1]-1)
            colB = colval_B[icolB]
            valB = nzval_B[icolB]
            for rowA in eachindex(rowptr_A[1:end-1])
                coliA = findfirst(col->col==rowB, colval_A[rowptr_A[rowA]:(rowptr_A[rowA+1]-1)])
                if(!isnothing(coliA))
                    coliA += rowptr_A[rowA]-1
                    valA = nzval_A[coliA]
                    C[(colB-1)*size(C,1) + rowA] = C[(colB-1)*size(C,1) + rowA] ⊕ op(valA, valB)
                end
            end
        end
    end

    return C
end

#= General parallel sparse product =#

#FIXME: Make sure all these functions are as warp-friendly as possible

function _spmdv_csc!(y::CuArray, op::Function, colptr_A::CuArray, rowval_A::CuArray, nzval_A::CuArray, x::CuArray) #FIXME: not sure if :: or <:

    y = fill!(y, zero(eltype(y)))

    kernel = @cuda launch=false _k_spmdv_csc!(y, op, colptr_A, rowval_A, nzval_A, x)

    threads, blocks = get_config(kernel, length(colptr_A)-1)

    CUDA.@sync kernel(y, op, colptr_A, rowval_A, nzval_A, x; threads=threads, blocks=blocks)

    return y
end

@inbounds function _k_spmdv_csc!(y, op, colptr_A, rowval_A, nzval_A, x)

    idx = get_idx()

    if idx <= length(colptr_A)-1
        for ptr in colptr_A[idx]:(colptr_A[idx+1]-1)
            CUDA.@atomic y[rowval_A[ptr]] = y[rowval_A[ptr]] ⊕ op(nzval_A[ptr], x[idx])
        end
    end
    nothing
end

function _spmdv_csr!(y::CuArray, op::Function, rowptr_A::CuArray, colval_A::CuArray, nzval_A::CuArray, x::CuArray)
    y = fill!(y, zero(eltype(y)))

    kernel = @cuda launch=false _k_spmdv_csr!(y, op, rowptr_A, colval_A, nzval_A, x)

    threads, blocks = get_config(kernel, length(rowptr_A)-1)

    CUDA.@sync kernel(y, op, rowptr_A, colval_A, nzval_A, x; threads=threads, blocks=blocks)

    return y
end

@inbounds function _k_spmdv_csr!(y, op, rowptr_A, colval_A, nzval_A, x)
    idx = get_idx()

    #the loop can probably be multi-threaded, tho it would lead to concurrent access anyway, so the time gain probably isn't very high
    if idx <= length(rowptr_A)-1
        for ptr in rowptr_A[idx]:(rowptr_A[idx+1]-1)
            CUDA.@atomic y[idx] = y[idx] ⊕ op(nzval_A[ptr], x[colval_A[ptr]])
        end
    end

    nothing
end

function _spmspv_csc!(y::CuArray, op::Function, colptr_A::CuArray, rowval_A::CuArray, nzval_A::CuArray, x_nzind::CuArray, x_nzval::CuArray)
    y = fill!(y, zero(eltype(y)))

    kernel = @cuda launch=false _k_spmspv_csc!(y, op, colptr_A, rowval_A, nzval_A, x_nzind, x_nzval)

    threads, blocks = get_config(kernel, length(x_nzind))

    CUDA.@sync kernel(y, op, colptr_A, rowval_A, nzval_A, x_nzind, x_nzval; threads=threads, blocks=blocks)

    return y
end

@inbounds function _k_spmspv_csc!(y, op, colptr_A, rowval_A, nzval_A, x_nzind, x_nzval)
    idx = get_idx()

    if idx <= length(x_nzind)
		xrow = x_nzind[idx]
        for ptr in colptr_A[xrow]:(colptr_A[xrow+1]-1)
            CUDA.@atomic y[rowval_A[ptr]] = y[rowval_A[ptr]] ⊕ op(nzval_A[ptr], x_nzval[idx])
        end
    end

    nothing
end

function _spmspv_csr!(y::CuArray, op::Function, rowptr_A::CuArray, colval_A::CuArray, nzval_A::CuArray, x_nzind::CuArray, x_nzval::CuArray)
    y = fill!(y, zero(eltype(y)))

    kernel = @cuda launch=false _k_spmspv_csr!(y, op, rowptr_A, colval_A, nzval_A, x_nzind, x_nzval)

    threads, blocks = get_config(kernel, length(rowptr_A)-1)

    CUDA.@sync kernel(y, op, rowptr_A, colval_A, nzval_A, x_nzind, x_nzval; threads=threads, blocks=blocks)

    return y
end

@inbounds function _k_spmspv_csr!(y, op, rowptr_A, colval_A, nzval_A, x_nzind, x_nzval)
    idx = get_idx()

    if idx <= length(rowptr_A)-1
        for ptr in rowptr_A[idx]:(rowptr_A[idx+1]-1)
            i = colval_A[ptr]

            j = findfirst(a -> a==i, x_nzind)
            if(!isnothing(j))
                CUDA.@atomic y[idx] = y[idx] ⊕ op(nzval_A[ptr], x_nzval[j])
            end
        end
    end

    nothing
end

#FIXME: I unplugged my brain for the last 2, it's probably not very efficient, not very warp friendly, etc...

function _spmspm_csc!(C::CuArray, op::Function, colptr_A::CuArray, rowval_A::CuArray, nzval_A::CuArray, colptr_B::CuArray, rowval_B::CuArray, nzval_B::CuArray)
    C = fill!(C, zero(eltype(C)))

    kernel = @cuda launch=false _k_spmspm_csc!(C, op, colptr_A, rowval_A, nzval_A, colptr_B, rowval_B, nzval_B)

    threads, blocks = get_config(kernel, length(colptr_B)-1)

    CUDA.@sync kernel(C, op, colptr_A, rowval_A, nzval_A, colptr_B, rowval_B, nzval_B; threads=threads, blocks=blocks)


    #set_kernel(_k_spmspm_csc!, length(colptr_A)-1 , C, op, colptr_A, rowval_A, nzval_A, colptr_B, rowval_B, nzval_B)

    return C
end

#TODO: Further testing required

@inbounds function _k_spmspm_csc!(C, op, colptr_A, rowval_A, nzval_A, colptr_B, rowval_B, nzval_B)
    idx = get_idx()

    if idx <= length(colptr_B)-1
        ptrrange = colptr_B[idx]:(colptr_B[idx+1]-1)
		_spmspv_csc!(view(C, :, idx), op, colptr_A, rowval_A, nzval_A, view(rowval_B, ptrrange), view(nzval_B, ptrrange))
    end

	nothing
end


function _spmspm_csr!(C::CuArray, op::Function, rowptr_A::CuArray, colval_A::CuArray, nzval_A::CuArray, rowptr_B::CuArray, colval_B::CuArray, nzval_B::CuArray)
    C = fill!(C, zero(eltype(C)))

    kernel = @cuda launch=false _k_spmspm_csr!(C, op, rowptr_A, colval_A, nzval_A, rowptr_B, colval_B, nzval_B)

    threads, blocks = get_config(kernel, length(rowptr_B)-1)

    CUDA.@sync kernel(C, op, rowptr_A, colval_A, nzval_A, rowptr_B, colval_B, nzval_B; threads=threads, blocks=blocks)

    #set_kernel(_k_spmspm_csr!, length(rowptr_B)-1, C, op, rowptr_A, colval_A, nzval_A, rowptr_B, colval_B, nzval_B)

    return C
end


#Could probably be improved

@inbounds function _k_spmspm_csr!(C, op, rowptr_A, colval_A, nzval_A, rowptr_B, colval_B, nzval_B)
	idx = get_idx()

    if idx <= length(rowptr_B)-1
        for ptrB in rowptr_B[idx]:(rowptr_B[idx+1]-1)
            colB = colval_B[ptrB]
            valB = nzval_B[ptrB]

            for rowA in 1:(length(rowptr_A)-1)
                ptrA = rowptr_A[rowA]:(rowptr_A[rowA+1]-1)
                coliA = findfirst(col->col==idx, view(colval_A, ptrA))

                if(!isnothing(coliA))
                    coliA += rowptr_A[rowA]-1
                    valA = nzval_A[coliA]
                    C[(colB-1)*size(C,1) + rowA] = C[(ptrB-1)*size(C,1) + rowA] ⊕ op(valA, valB)
                end
            end
        end
    end
	nothing
end
