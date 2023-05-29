### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ e2560be8-f3f9-11ed-0754-e3e1572c751d
begin
	using Pkg
	#Pkg.activate(mktempdir())
	
	# Register the FAST registry
	Pkg.Registry.add(
		RegistrySpec(url="https://gitlab.lisn.upsaclay.fr/fast/registry")
	)

	using Revise
	Pkg.develop(path="..")

	using FiniteStateAutomata

	Pkg.add("Semirings")
	using Semirings 

	using LinearAlgebra
end

# ╔═╡ 2b8f46d1-08a5-479a-8261-7dc6f662563a
K = LogSemiring{Float64}

# ╔═╡ b0f12741-f4b1-48f5-bf75-756b55d99975
md"""## Closure"""

# ╔═╡ fc2aedf5-3050-4120-a138-efd10be3994d
T = TransitionMatrix(
	K, 
	3,
	[
		(1, 1, K(2)), 
		(1, 2, K(3)), 
		(2, 3, K(4)), 
		(3, 3, K(4))
	],
	[(1, 1, K(2))],
	[(1, 2, K(2))],
	[(1, 2, K(3)), (2, 3, K(3))]
);

# ╔═╡ 5ea5e692-1e92-4799-b355-f19c35645341
A = FST(
	K,
	[(1, K(5))],
	T,
	[(2, K(5)), (3, K(5))],
	["a" => "p", "b" => "q", "c" => "r"],
	# infactors = [(1, 1, K(2))],
	# outfactors = [(1, 2, K(3)), (1, 3, K(3))]
	# factors = [(1, 2, K(2)]
	# [2, []]
)

# ╔═╡ 2e4cee38-1cf9-499b-8ec8-5f92e0286b4e
closure(A) 

# ╔═╡ de8cc675-45d6-4d7b-84b1-b93dcfdccebf
A

# ╔═╡ 04765157-f208-4e04-850d-3e83d7fa28f3
closure(A) |> T

# ╔═╡ 61f9e6fc-96a7-4ba1-99e6-89b1303552f6
A₁ = FST(
	K,
	[(0, 1, one(K)), (0, 2, one(K)), (1, 3, one(K)), (1, 2, one(K)), (2, 3, one(K))],
	[(3, one(K))],
	["a" => "x", "b" => "y", "c" => "z"]
)

# ╔═╡ 797af1e7-40d4-4fc6-aebc-fc8daffd232e
hcat(α(A₁)) * vcat(transpose(ω(A₁)))

# ╔═╡ b774dd87-a29b-4fbe-966f-0b797199a735
A₁  |> determinize

# ╔═╡ e20b21d2-5700-4540-9d96-f1634c198587
union(A₁, A₁)  

# ╔═╡ f8bdd2a9-7034-4deb-a52b-f7402d05c7dc
transpose(α(A₁)) * T(A₁)

# ╔═╡ 7da1862a-3e7c-4d8f-ab85-7c115871b67c
closure(A₁; plus = false)

# ╔═╡ 17734edd-7c63-4f3f-8800-942d928629f4
A₁ |> renorm

# ╔═╡ 77e1dc3b-ca61-48a1-8324-6a974a69b950
ρ(A₁) * α(A₁)

# ╔═╡ 0edfc9a3-31d5-47b3-817a-28ef40c69f9a
reverse(A₁) 

# ╔═╡ 27b3d329-c2c1-4006-a7fa-4ec37b61c677
sum(eachcol(T(A₁)))

# ╔═╡ 32ccbf5c-f3bf-43d0-a89f-8b27a18a44e7
W(A₁)

# ╔═╡ e3d2e145-3349-40e7-b24e-5133fbe2c0a7
A₂ = FST(
	K,
	[(0, 1, one(K)), (1, 2, one(K)), (2, 1, one(K)), (2, 3, one(K))],
	[(3, one(K))],
	["a" => "u", "b" => "v", "c" => "w"]
)

# ╔═╡ 80db008d-1be1-4411-b1ad-5e0757b5750b
A₁ ∪ A₂

# ╔═╡ f676d57d-31dc-4861-8dcd-e12eeeba48d4
cat(A₁, A₂) |> reverse 

# ╔═╡ eb907a8c-9233-4615-9e73-3297f4f48535
Iterators.map(x -> 2x, [1, 2, 3]) |> typeof

# ╔═╡ Cell order:
# ╠═e2560be8-f3f9-11ed-0754-e3e1572c751d
# ╠═2b8f46d1-08a5-479a-8261-7dc6f662563a
# ╟─b0f12741-f4b1-48f5-bf75-756b55d99975
# ╠═fc2aedf5-3050-4120-a138-efd10be3994d
# ╠═5ea5e692-1e92-4799-b355-f19c35645341
# ╠═2e4cee38-1cf9-499b-8ec8-5f92e0286b4e
# ╠═de8cc675-45d6-4d7b-84b1-b93dcfdccebf
# ╠═04765157-f208-4e04-850d-3e83d7fa28f3
# ╠═797af1e7-40d4-4fc6-aebc-fc8daffd232e
# ╠═61f9e6fc-96a7-4ba1-99e6-89b1303552f6
# ╠═80db008d-1be1-4411-b1ad-5e0757b5750b
# ╠═b774dd87-a29b-4fbe-966f-0b797199a735
# ╠═e20b21d2-5700-4540-9d96-f1634c198587
# ╠═f8bdd2a9-7034-4deb-a52b-f7402d05c7dc
# ╠═f676d57d-31dc-4861-8dcd-e12eeeba48d4
# ╠═7da1862a-3e7c-4d8f-ab85-7c115871b67c
# ╠═17734edd-7c63-4f3f-8800-942d928629f4
# ╠═77e1dc3b-ca61-48a1-8324-6a974a69b950
# ╠═0edfc9a3-31d5-47b3-817a-28ef40c69f9a
# ╠═27b3d329-c2c1-4006-a7fa-4ec37b61c677
# ╠═32ccbf5c-f3bf-43d0-a89f-8b27a18a44e7
# ╠═e3d2e145-3349-40e7-b24e-5133fbe2c0a7
# ╠═eb907a8c-9233-4615-9e73-3297f4f48535
