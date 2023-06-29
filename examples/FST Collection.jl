### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ a183e6f2-157e-11ee-0d03-ad978a190e29
begin
	using Pkg
	Pkg.develop(path="../")

	using Revise
	using FiniteStateAutomata
end

# ╔═╡ da6fdaa9-675d-49b8-80e7-a5aaae988b10
symsdir = "../examples/symboltables/"

# ╔═╡ 7931181b-f0dd-4d9c-b5cf-a7f1a073810f
fstdir = "../examples/fsts/"

# ╔═╡ 57178d12-b53c-4207-a66f-21d3d0c1288b
symtables = Dict(
	:ascii => open(loadsymbols, joinpath(symsdir, "ascii.syms")),
	:latin => open(loadsymbols, joinpath(symsdir, "latin.syms")),
	:emoticons => open(loadsymbols, joinpath(symsdir, "emoticons.syms")),
	:cyrillic => open(loadsymbols, joinpath(symsdir, "cyrillic.syms"))
)

# ╔═╡ d0b1040d-dbd0-4073-b9e6-6548762fd955
open(compile, joinpath(fstdir, "fst_ex1.txt"))

# ╔═╡ Cell order:
# ╠═a183e6f2-157e-11ee-0d03-ad978a190e29
# ╠═da6fdaa9-675d-49b8-80e7-a5aaae988b10
# ╠═7931181b-f0dd-4d9c-b5cf-a7f1a073810f
# ╠═57178d12-b53c-4207-a66f-21d3d0c1288b
# ╠═d0b1040d-dbd0-4073-b9e6-6548762fd955
