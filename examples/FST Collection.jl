### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ a183e6f2-157e-11ee-0d03-ad978a190e29
begin
	using Pkg
	Pkg.activate("./env")
	Pkg.develop(path="../")

	using Revise
	using FiniteStateAutomata
	using NPZ
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
S = ProbSemiring{Float32}

# ╔═╡ f6598817-db39-43d4-8614-86bf47e30f9f
vectorfst = VectorFST(
	[
		[(2, 1, 1, S(.5)), (3, 2, 2, S(1.5))],
		[(3, 3, 3, S(2.5))],
		Arc{S}[]
	],
	1,
	S[zero(S), zero(S), S(3.5)]
)

# ╔═╡ 2fc5a5c3-7d40-4b6b-a7da-a338608ffd25
t = reshape([S(0.1) for i in 1:36], (2, 2, 3, 3))

# ╔═╡ 19e08213-151a-4945-a989-53e0b369d3e5


# ╔═╡ 6a70a757-f606-41f4-941d-49c7e2b814d8
size(t)

# ╔═╡ 37d55fc1-dcdc-44d0-9a70-b3b4d7d796c5
#Transducer for 2 states 3 labels
tensorFST = TensorFST(
	t,
	[S(3.5), zero(S)],
	[zero(S), S(2.0)]

)

# ╔═╡ 8548092c-7902-4c38-b007-6799b6a346f8
setinitstate!(vectorfst, 3)

# ╔═╡ 3739bac6-b90e-444a-aa6f-03dd7f1e765b
draw(
	tensorFST; 
	isymbols=symtables[:latin], 
	osymbols=symtables[:cyrillic]
) 

# ╔═╡ 8b543c88-4136-4c77-935f-0e617f273a65
md"""
## Example 3: NN outputs as DenseFST
"""

# ╔═╡ f67183da-92b1-421b-834e-bd4f3dd815d4
logits = NPZ.npzread("./assets/libri_examples/2830-3980-0002/logits.npy")

# ╔═╡ a8a3c4ec-ca78-4b0e-ad8f-338b48f5c956
begin
	label_mapping = open("./assets/libri_examples/ctc_map.txt") do f
		readlines(f)
	end
	label_mapping = Dict(k => v for (k,v) in enumerate(label_mapping))
end

# ╔═╡ 8dbc2622-22b6-49df-8e7e-3c7e2e10e01b
begin
	# lets represent the logits as label-dependent alignment lattice
	local Q, L = size(logits)
	Q += 1 # +1 for initial/final state
	local K = LogSemiring{Float16, ℯ}
	
	local W = K.(logits)
	local α = zeros(K, Q)
	α[1] = one(K)
	local ω = zeros(K, Q)
	ω[end] = one(K)
	local M = zeros(K, Q, Q, L, L)

	for q in 1:(Q-1)
		for l in 1:L
			M[q, q+1, l, l] = W[q, l]
 		end
	end
	nn_fst = TensorFST(M, α, ω)
	draw(nn_fst; isymbols=label_mapping, osymbols=label_mapping)
end

# ╔═╡ Cell order:
# ╠═a183e6f2-157e-11ee-0d03-ad978a190e29
# ╠═da6fdaa9-675d-49b8-80e7-a5aaae988b10
# ╠═7931181b-f0dd-4d9c-b5cf-a7f1a073810f
# ╠═57178d12-b53c-4207-a66f-21d3d0c1288b
# ╠═d0b1040d-dbd0-4073-b9e6-6548762fd955
# ╠═5517d174-48cf-4c33-9646-e52c57c1ac4e
# ╠═f6598817-db39-43d4-8614-86bf47e30f9f
# ╠═2fc5a5c3-7d40-4b6b-a7da-a338608ffd25
# ╠═6a70a757-f606-41f4-941d-49c7e2b814d8
# ╠═37d55fc1-dcdc-44d0-9a70-b3b4d7d796c5
# ╠═8548092c-7902-4c38-b007-6799b6a346f8
# ╠═3739bac6-b90e-444a-aa6f-03dd7f1e765b
# ╟─8b543c88-4136-4c77-935f-0e617f273a65
# ╠═f67183da-92b1-421b-834e-bd4f3dd815d4
# ╠═a8a3c4ec-ca78-4b0e-ad8f-338b48f5c956
# ╠═8dbc2622-22b6-49df-8e7e-3c7e2e10e01b
