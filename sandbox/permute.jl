### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ a183e6f2-157e-11ee-0d03-ad978a190e29
begin
	using Pkg
	Pkg.activate("../")

	using Revise
	using FiniteStateAutomata
end

# ╔═╡ dcbb87c3-139b-465c-a1e0-04774a868411
S = LogSemiring{Float32,1}

# ╔═╡ da6fdaa9-675d-49b8-80e7-a5aaae988b10
symsdir = "../examples/symboltables/"

# ╔═╡ 49fead47-9d06-4919-a675-ea1591aa83eb
syms = open(loadsymbols, "../examples/symboltables/ascii.syms")

# ╔═╡ 96d478e3-f21f-4fa8-ad76-3f88424d3a52
fst = convert(
	TensorFST{S,Array{S,4}},
	open("../examples/fsts/fst_ex1.txt") do f
		compile(f; semiring = S)
	end
)

# ╔═╡ b31e8f95-1218-4a20-873a-c8812b306270


# ╔═╡ Cell order:
# ╠═a183e6f2-157e-11ee-0d03-ad978a190e29
# ╠═dcbb87c3-139b-465c-a1e0-04774a868411
# ╠═da6fdaa9-675d-49b8-80e7-a5aaae988b10
# ╠═49fead47-9d06-4919-a675-ea1591aa83eb
# ╠═96d478e3-f21f-4fa8-ad76-3f88424d3a52
# ╠═b31e8f95-1218-4a20-873a-c8812b306270
