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
#open(compile, joinpath(fstdir, "fst_ex1.txt"))

# ╔═╡ 5517d174-48cf-4c33-9646-e52c57c1ac4e
S = LogSemiring{Float32,1}

# ╔═╡ f6598817-db39-43d4-8614-86bf47e30f9f
vectorfst = VectorFST(
	[
		[(2, 1 => 1, S(.5)), (3, 2 => 2, S(1.5))],
		[(3, 3 => 3, S(2.5))],
		Tuple{Int,Pair{Int,Int},S}[]
	],
	1,
	S[zero(S), zero(S), S(3.5)]
)

# ╔═╡ 8548092c-7902-4c38-b007-6799b6a346f8
setinitstate!(vectorfst, 3)

# ╔═╡ 3739bac6-b90e-444a-aa6f-03dd7f1e765b
draw(
	vectorfst; 
	isymbols=symtables[:latin], 
	osymbols=symtables[:cyrillic]
) 

# ╔═╡ 595a866c-635e-45bb-839e-fcb6ad6f396e
sum([length(collect(arcs(vectorfst, q))) for q in states(vectorfst)])

# ╔═╡ 45085224-629a-4596-8741-56d4c941e762
draw(
	deletestate!(vectorfst, 2); 
	isymbols=symtables[:latin], 
	osymbols=symtables[:cyrillic]
) 

# ╔═╡ Cell order:
# ╠═a183e6f2-157e-11ee-0d03-ad978a190e29
# ╠═da6fdaa9-675d-49b8-80e7-a5aaae988b10
# ╠═7931181b-f0dd-4d9c-b5cf-a7f1a073810f
# ╠═57178d12-b53c-4207-a66f-21d3d0c1288b
# ╠═d0b1040d-dbd0-4073-b9e6-6548762fd955
# ╠═5517d174-48cf-4c33-9646-e52c57c1ac4e
# ╠═f6598817-db39-43d4-8614-86bf47e30f9f
# ╠═8548092c-7902-4c38-b007-6799b6a346f8
# ╠═3739bac6-b90e-444a-aa6f-03dd7f1e765b
# ╠═595a866c-635e-45bb-839e-fcb6ad6f396e
# ╠═45085224-629a-4596-8741-56d4c941e762
