# SPDX-License-Identifier: CECILL-2.1

#=

Helper functions for operations on GPU.

=#

# Helper function to get the number of blocks and threads for a given kernel iterating on an array of length n
function get_config(k::CUDA.HostKernel, n)
    threads = min(n, launch_configuration(k.fun).threads)
    blocks = cld(n, threads)
    return threads, blocks
end

# Helper function to get the 1D index of a thread, amongst multiple blocks
get_idx() = threadIdx().x + (blockIdx().x - 1)*blockDim().x

function set_kernel(f::Function, lim::Integer, args::Vararg)
    kernel = @cuda launch=false f(args)

    threads, blocks = get_config(kernel, lim)

    CUDA.@sync kernel(args; threads=threads, blocks=blocks)
end

# Atomic defined for Semirings
@inline function CUDA.atomic_cas!(ptr::Core.LLVMPtr{T,A}, cmp::T, new::T) where {T<:Semiring, A}
	IT = CUDA.inttype(typeof(val(cmp)))
	cmp_i = reinterpret(IT, val(cmp))
	new_i = reinterpret(IT, val(new))
	old_i = CUDA.atomic_cas!(reinterpret(CUDA.LLVMPtr{IT,A}, ptr), cmp_i, new_i)
	return T(reinterpret(typeof(val(cmp)), old_i))
end

# Reduction over a warp
function warp_reduce(op::Function, x::T) where T <:Semiring
    offset = warpsize() ÷ 2
    while offset > 0
        x = op(x, T(CUDA.shfl_down_sync(CUDA.FULL_MASK, val(x), offset)))
        offset ÷= 2
    end
    x
end

function warp_reduce(op::Function, x::T) where T
    offset = warpsize() ÷ 2
    while offset > 0
        x = op(x, CUDA.shfl_down_sync(CUDA.FULL_MASK, x, offset))
        offset ÷= 2
    end
    x
end

# Reduction over all blocks
function all_reduce!(acc::CUDA.CuDeviceVector, op::Function, x::T, lane::Int) where T
    x = warp_reduce(op, x)

    if lane == 1
        @inbounds CUDA.@atomic acc[1] = op(acc[1], x)
    end
end

