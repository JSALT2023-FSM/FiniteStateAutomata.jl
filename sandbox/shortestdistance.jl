### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ f939691a-14e0-11ee-0cd4-83adb8b7a52c
begin
	using Pkg
	Pkg.develop(path="..")

	using Revise
	using FiniteStateAutomata
end

# ╔═╡ b12ab39b-7d64-43ab-b8c1-ba8c14e081c8
S = ProbSemiring{Float32}

# ╔═╡ ad860423-17fc-4cbf-9838-00685685ea7c
A = SparseFST(
	SparseMatrices(
		sparse([1, 2], [2, 3], one(S), 4, 4),
		sparse([1], [3], one(S), 4, 4),
		sparse([2, 3], [4, 4], one(S), 4, 4),
	),
	sparsevec([1], one(S), 4),
	sparsevec([3, 4], one(S), 4),
	[1, 2, 3]
)

# ╔═╡ f205a6bf-0aff-43f8-9850-1fc2f9c2eb16
T = sum(M(A)) 

# ╔═╡ 0dc9667c-4ef8-485e-bab4-a88b73147d8c
Tstar = powerseries(T);

# ╔═╡ f94860b0-99a6-4c30-92d9-a33121fdac65
draw(A; symbols=Dict(1 => "a", 2 => "b", 3 => "c"))

# ╔═╡ 0d1fa1bf-7c40-4060-96f7-13c585faeb51
transpose(α(A)) * Tstar

# ╔═╡ 5a407799-f256-4db4-9cf0-09efc86b927d
x = Array(α(A))

# ╔═╡ 47551543-6803-4aaf-997f-cb55e15d301d
B = Array(T)

# ╔═╡ aea8480c-e1e0-417b-9c32-0fb55e1db0f1
transpose(x) * B * B * B * B 

# ╔═╡ Cell order:
# ╠═f939691a-14e0-11ee-0cd4-83adb8b7a52c
# ╠═b12ab39b-7d64-43ab-b8c1-ba8c14e081c8
# ╠═ad860423-17fc-4cbf-9838-00685685ea7c
# ╠═f205a6bf-0aff-43f8-9850-1fc2f9c2eb16
# ╠═0dc9667c-4ef8-485e-bab4-a88b73147d8c
# ╠═f94860b0-99a6-4c30-92d9-a33121fdac65
# ╠═0d1fa1bf-7c40-4060-96f7-13c585faeb51
# ╠═5a407799-f256-4db4-9cf0-09efc86b927d
# ╠═47551543-6803-4aaf-997f-cb55e15d301d
# ╠═aea8480c-e1e0-417b-9c32-0fb55e1db0f1
