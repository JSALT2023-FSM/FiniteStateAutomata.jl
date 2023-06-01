### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ e2560be8-f3f9-11ed-0754-e3e1572c751d
begin
	using Pkg
	Pkg.activate(mktempdir())
	
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
	using SparseArrays
end

# ╔═╡ 2b8f46d1-08a5-479a-8261-7dc6f662563a
K = LogSemiring{Float64}

# ╔═╡ b0f12741-f4b1-48f5-bf75-756b55d99975
md"""## Closure"""

# ╔═╡ 3b996f92-639b-4f08-baa8-040db7190657
WFST(
	semiring = K,
	initweights = [(1, K(5))],
	arcs = [(1, 3, K(2)), (1, 2, K(3)), (2, 3, K(4))],
	finalweights = [(2, K(5)), (3, K(5))],
	ϵweight = one(K),
	statelabels = ["a" => "x", "b" => "y", "c" => "z"],
)

# ╔═╡ 4697649b-2925-4e90-b72e-08f42b975401
WFST(
	semiring = K,
	initweights = [(1, K(5))],
	arcs = [(1, 3, K(2)), (1, 2, K(3)), (2, 3, K(4))],
	finalweights = [(2, K(5)), (3, K(5))],
	statelabels = ["a", "b", "c"],
)

# ╔═╡ 97d94beb-eff7-44d0-b29b-966caab19470


# ╔═╡ 5ea5e692-1e92-4799-b355-f19c35645341
A = WFST(
	semiring = K,
	initweights = [(1, K(5))],
	arcs = [(1, 3, K(2)), (1, 2, K(3)), (2, 3, K(4))],
	finalweights = [(2, K(5)), (3, K(5))],
	statelabels = ["a" => "x", "b" => "y", "c" => "z"],
	infactors = [(1, 1, K(2))],
	outfactors = [(1, 2, K(3)), (2, 3, K(3))],
	factors = [(1, 2, K(2))]
)

# ╔═╡ 996ca0c1-ce6d-4f12-b52c-f0745f274b2c
Π₂(A)

# ╔═╡ 04cbae24-d6c3-462f-9665-fc4edd72e3d1
Π₂(A) isa FiniteStateAutomata.AbstractWFST{<:Semiring,<:String}

# ╔═╡ 2bb54fe3-e5f4-4dca-943f-235865a45555
B = FST(
	semiring = K,
	initweights = [(1, K(5))],
	arcs = [
		(1, 3, K(2)), 
		(1, 2, K(3)), 
		(2, 3, K(4))
	],
	finalweights = [(2, K(5)), (3, K(5))],
	statelabels = ["a" => "x", "b" => "y", "c" => "z"],
)

# ╔═╡ e4fadd4b-6e11-44c7-b4ed-b7add8c9071a
B isa Transducer

# ╔═╡ b2e09b15-8bc7-4503-a37f-fceeae5e14bc
FiniteStateAutomata.MatrixPowerSum(T(A)) * ω(A)

# ╔═╡ 2e4cee38-1cf9-499b-8ec8-5f92e0286b4e
closure(A; plus = true) 

# ╔═╡ de8cc675-45d6-4d7b-84b1-b93dcfdccebf
closure(B; plus = true)

# ╔═╡ 797af1e7-40d4-4fc6-aebc-fc8daffd232e
A ∪ A

# ╔═╡ 256a5845-e828-4490-b495-f5a7753d2375
B ∪ B ∪ B

# ╔═╡ 7dfd3464-78ac-4f5a-ab1c-ccfe009d9cd6
cat(A, A)

# ╔═╡ ce942c70-2f46-4b9c-9264-c83f1225dfc9
cat(B, B)

# ╔═╡ bd07d3fa-35d3-4f21-a136-8cd11aea6db7
Sa = sprandn(3, 3, 0.2)

# ╔═╡ 2bc05444-e9dd-40c4-9994-aa35de93086f
Sb = sprandn(3, 3, 0.2)

# ╔═╡ d09c6abe-32a9-49af-b39d-ba2323f484d2
x = sprandn(9, 0.2)

# ╔═╡ c696c786-3ddb-4c24-a44c-c0de1d4800c0
kron(Sa, Sb) * x

# ╔═╡ 0e589389-26e2-497a-ae70-9bd85ed6881b
vec(Sb * reshape(x, 3, 3) * transpose(Sa))

# ╔═╡ 61f9e6fc-96a7-4ba1-99e6-89b1303552f6
A₁ = FST(
	K,
	[(0, 1, one(K)), (0, 2, one(K)), (1, 3, one(K)), (1, 2, one(K)), (2, 3, one(K))],
	[(3, one(K))],
	["a" => "x", "b" => "y", "c" => "z"]
)

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

# ╔═╡ f676d57d-31dc-4861-8dcd-e12eeeba48d4
cat(A₁, A₂) |> reverse 

# ╔═╡ eb907a8c-9233-4615-9e73-3297f4f48535
Iterators.map(x -> 2x, [1, 2, 3]) |> typeof

# ╔═╡ Cell order:
# ╠═e2560be8-f3f9-11ed-0754-e3e1572c751d
# ╠═2b8f46d1-08a5-479a-8261-7dc6f662563a
# ╟─b0f12741-f4b1-48f5-bf75-756b55d99975
# ╠═3b996f92-639b-4f08-baa8-040db7190657
# ╠═4697649b-2925-4e90-b72e-08f42b975401
# ╠═97d94beb-eff7-44d0-b29b-966caab19470
# ╠═5ea5e692-1e92-4799-b355-f19c35645341
# ╠═996ca0c1-ce6d-4f12-b52c-f0745f274b2c
# ╠═04cbae24-d6c3-462f-9665-fc4edd72e3d1
# ╠═2bb54fe3-e5f4-4dca-943f-235865a45555
# ╠═e4fadd4b-6e11-44c7-b4ed-b7add8c9071a
# ╠═b2e09b15-8bc7-4503-a37f-fceeae5e14bc
# ╠═2e4cee38-1cf9-499b-8ec8-5f92e0286b4e
# ╠═de8cc675-45d6-4d7b-84b1-b93dcfdccebf
# ╠═797af1e7-40d4-4fc6-aebc-fc8daffd232e
# ╠═256a5845-e828-4490-b495-f5a7753d2375
# ╠═7dfd3464-78ac-4f5a-ab1c-ccfe009d9cd6
# ╠═ce942c70-2f46-4b9c-9264-c83f1225dfc9
# ╠═bd07d3fa-35d3-4f21-a136-8cd11aea6db7
# ╠═2bc05444-e9dd-40c4-9994-aa35de93086f
# ╠═d09c6abe-32a9-49af-b39d-ba2323f484d2
# ╠═c696c786-3ddb-4c24-a44c-c0de1d4800c0
# ╠═0e589389-26e2-497a-ae70-9bd85ed6881b
# ╠═61f9e6fc-96a7-4ba1-99e6-89b1303552f6
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
