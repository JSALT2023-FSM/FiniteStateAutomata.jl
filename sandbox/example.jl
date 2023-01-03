### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ e783cfcc-8b82-11ed-147f-d18598bc60ff
begin
	using Pkg
	Pkg.activate("../")

	using Revise
	using FSALinalg
	using LinearAlgebra
	using SparseArrays
end

# ╔═╡ 87e1cf8e-da52-4c85-ab6a-bb352a84d779
K = Float32

# ╔═╡ 8602a2b5-16ca-4cec-8041-c6362aa66a48
FSA(
	sparsevec([1, 2], K[0.5, 1.5], 3),
	sparse([1, 2], [3, 3], K[2.5, 2.5], 3, 3),
	sparsevec([3], K[3.5], 3),
	["a", "b", "c"]
)

# ╔═╡ Cell order:
# ╠═e783cfcc-8b82-11ed-147f-d18598bc60ff
# ╠═87e1cf8e-da52-4c85-ab6a-bb352a84d779
# ╠═8602a2b5-16ca-4cec-8041-c6362aa66a48
