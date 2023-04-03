### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 510f7d8a-ce11-11ed-2832-cdc20fa6f700
begin
	using Pkg; Pkg.activate("..")
	using Revise
	using FiniteStateAutomata
	using Semirings
	using Random
	using SparseArrays
	using LinearAlgebra
end

# ╔═╡ afc66d19-c7c3-4318-8fa0-f5fc640ca271
K = LogSemiring{Float32}

# ╔═╡ 2c44097b-35f7-4d50-9e27-a703c653ceb0
Σ = Set(["a", "b", "c"])

# ╔═╡ fe9d52c8-c152-4961-847a-ac9a1554e07b
H = K.(randn(length(Σ), 4))

# ╔═╡ cf73fdb5-b276-4f4f-b174-cf4db4f4d216
N = sparse([], [], K[], length(Σ), length(Σ))

# ╔═╡ ab6a60d4-8c71-4050-9d3d-37e49128da25
e = ones(K, length(Σ))

# ╔═╡ ba26ff62-d337-418b-8992-8230484b8cab
Q = size(H, 1) * (size(H, 2) - 1)

# ╔═╡ 313206a0-b8cb-4766-9ac1-fe4ddd29e8d2
L = FSA(
	sparsevec(1:length(Σ), H[:, 1], Q),
	sparse([
		N sparse(e * H[:, 2]') N;
		N N sparse(e * H[:, 3]');
		N N N
	]),
	sparsevec(Q - length(Σ) + 1: Q, H[:, end], Q),
	zero(K),
	repeat(["a", "b", "c"], 3)
)

# ╔═╡ d123178c-561d-4780-bde8-6edcc21f9ae9
N

# ╔═╡ 66659283-01c5-4883-854a-5a03012fb56e
C = FSA(
	sparsevec([1], one(K), 4),
	sparse([1, 1, 1, 2, 3, 4], [1, 2, 3, 4, 4, 4], one(K), 4, 4),
	sparsevec([4], one(K), 4),
	zero(K),
	["a", "b", "a", "c"]
)

# ╔═╡ 770b840b-c9cd-470b-bb7c-861a33e4fb14
D = FSA(
	sparsevec([1], one(K), 3),
	sparse([1, 1, 2, 3], [1, 2, 3, 1], one(K), 3, 3),
	sparsevec([3], one(K), 3),
	zero(K),
	["a", "b", "c"]
)

# ╔═╡ 7820e2e6-55df-49ee-8fde-ffb4c69269b5
C ∩ D

# ╔═╡ b37d3f48-0b23-4700-bb45-050b08209669
@time L ∩ C

# ╔═╡ 183f8031-b384-408b-8cd6-c84f23693dfc
# ╠═╡ disabled = true
#=╠═╡
M = sparse([1,1,2,3,4,4,5,6,7,7,8,9],[1,3,2,4,1,3,2,4,1,3,2,4], one(K), nstates(L), nstates(C))
  ╠═╡ =#

# ╔═╡ b605383a-4080-4916-857c-2528f6ad172c
M' * L.ω

# ╔═╡ e91b7104-7149-4770-8abe-fb0a968bd197
L

# ╔═╡ 49fb8a23-bae8-4c66-8490-01d7e28f25f9
begin
	e3 = ones(K, nstates(C))
	N3 = spzeros(K, nstates(C), nstates(C))
	M3 = sparse([1,1,2,3],[1,3,2,4], one(K), size(H,1), nstates(C))
	T3 = [
		N3 sparse(e3 * (M3'*H[:, 2])') .* C.T N3;
		N3 N3 sparse(e3 * (M3'*H[:, 2])') .* C.T;
		N3 N3 N3;
	]
	T3
end

# ╔═╡ 29a16022-6ebc-4e74-8245-e7293fcb160a
A3=FSA(
	sparsevec(1:nstates(C), C.α .* (M3' * H[:, 1]), nstates(C) * (size(H, 2)- 1)),
	T3,
	sparsevec( (size(H, 2) - 2) * nstates(C) + 1 : (size(H, 2) - 1) * nstates(C), C.ω .* (M3' * H[:,end]), nstates(C) * (size(H, 2)- 1)),
	C.ρ * L.ρ,
	repeat(C.λ, size(H, 2) -1)
)

# ╔═╡ ada373ae-0e4e-42ca-80a5-98be29aca610
G = DenseFSA(H, sort(collect(Σ)), zero(K))

# ╔═╡ 049d12b4-0107-4ecc-89e4-45d3955bdd0b
@time G ∩ C

# ╔═╡ 54888ef8-6ee9-446d-9b78-364fe35e52b6
G, L

# ╔═╡ 40ba3c52-84e0-4d70-8d8e-840ec2225baf
G ∩ C

# ╔═╡ 91b52b46-b8d8-4d14-9f9c-732de9d906ff
G ∩ C |> connect

# ╔═╡ Cell order:
# ╠═510f7d8a-ce11-11ed-2832-cdc20fa6f700
# ╠═afc66d19-c7c3-4318-8fa0-f5fc640ca271
# ╠═2c44097b-35f7-4d50-9e27-a703c653ceb0
# ╠═fe9d52c8-c152-4961-847a-ac9a1554e07b
# ╠═cf73fdb5-b276-4f4f-b174-cf4db4f4d216
# ╠═ab6a60d4-8c71-4050-9d3d-37e49128da25
# ╠═ba26ff62-d337-418b-8992-8230484b8cab
# ╠═313206a0-b8cb-4766-9ac1-fe4ddd29e8d2
# ╠═d123178c-561d-4780-bde8-6edcc21f9ae9
# ╠═66659283-01c5-4883-854a-5a03012fb56e
# ╠═770b840b-c9cd-470b-bb7c-861a33e4fb14
# ╠═7820e2e6-55df-49ee-8fde-ffb4c69269b5
# ╠═b37d3f48-0b23-4700-bb45-050b08209669
# ╠═049d12b4-0107-4ecc-89e4-45d3955bdd0b
# ╠═54888ef8-6ee9-446d-9b78-364fe35e52b6
# ╠═183f8031-b384-408b-8cd6-c84f23693dfc
# ╠═b605383a-4080-4916-857c-2528f6ad172c
# ╠═e91b7104-7149-4770-8abe-fb0a968bd197
# ╠═49fb8a23-bae8-4c66-8490-01d7e28f25f9
# ╠═29a16022-6ebc-4e74-8245-e7293fcb160a
# ╠═ada373ae-0e4e-42ca-80a5-98be29aca610
# ╠═40ba3c52-84e0-4d70-8d8e-840ec2225baf
# ╠═91b52b46-b8d8-4d14-9f9c-732de9d906ff
