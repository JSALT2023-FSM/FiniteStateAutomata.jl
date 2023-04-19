### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ 00f12bce-d9dc-11ed-2b99-edb4fd8c8b94
begin
	using Pkg
	Pkg.activate("../")

	using Revise
	using FiniteStateAutomata
	using SparseArrays
	using Semirings
	using Plots
	using Zygote

	using PlutoUI
end

# ╔═╡ 08803e3f-cf88-4d91-962f-220709e51019
md"""
# Taking derivatives with respect to an FST

This notebook is a minimal example of how to implement a FST-based loss function
and backpropagate through it.
"""

# ╔═╡ 1e8e64d2-38be-4474-af8e-077426bbdc92
md"""
First step is to define the semiring we are going to work on. Here we use the log-semiring.
"""

# ╔═╡ 86528707-e1b5-4cbc-bc06-270ead357ada
K = LogSemiring{Float32}

# ╔═╡ 2bde23e4-010f-4102-97b0-5500e1e85fb0
md"""
In our example we have 2 input FST `A1` and `A2`. In practice their weights would be given by a neural network.
"""

# ╔═╡ ccb356bb-9592-4f04-a4fe-79d7d2e439f4
A1 = FST(
	sparsevec([1, 2], one(K), 4),
	sparse([1, 2, 2], [3, 3, 4], one(K), 4, 4),
	sparsevec([3, 4], one(K), 4),
	zero(K),
	["a", "b", "c", "d"]
)

# ╔═╡ 29fb4884-84b5-4e54-a339-50703c1bb89f
A2 = FST(
	sparsevec([1, 2], one(K), 4),
	sparse([1, 1, 1, 2, 4], [2, 3, 4, 3, 2], one(K), 4, 4),
	sparsevec([3, 4], one(K), 4),
	zero(K),
	["a", "b", "d", "c"]
)

# ╔═╡ 70b8e8f7-e6c1-4dbc-b165-6eb301395852
md"""
We also have some "soft supervision", i.e. the supervision is not a unique sequence but rather a set of sequence each with a particular weighting represeneted as FST. In this example the supervision FST will be `B1` and `B2`.
"""

# ╔═╡ f7619ff9-7142-4964-b981-d6603bb3df07
B1 = FST(
	sparsevec([1], one(K), 3),
	sparse([1, 1, 2, 2], [2, 3, 1, 3], one(K), 3, 3),
	sparsevec([2, 3], one(K), 3),
	zero(K),
	["a", "b", "c"]
)

# ╔═╡ a09dd6a4-560c-403d-8278-819f16c8451a
B2 = FST(
	sparsevec([1, 2], one(K), 3),
	sparse([1, 1, 2, 2], [2, 3, 1, 3], one(K), 3, 3),
	sparsevec([2, 3], one(K), 3),
	zero(K),
	["a", "b", "c"]
)

# ╔═╡ 1c3a8f7e-71e4-4db7-8cc1-9a327d006780
C = FST(
	sparsevec(1:4, one(K), 4),
	sparse(ones(K, 4, 4) * K(log(1/4))),
	sparsevec(1:4, one(K), 4),
	zero(K),
	["a", "b", "c", "d"]
) 

# ╔═╡ 2ef2e1c1-6872-4d63-8d15-8b118cbad37f
begin
	X1 = FST(A1)
	X2 = FST(A2)
	local numiter = 200
	local τ = 1e-1

	F = zeros(numiter)
	for i in 1:numiter
		(F[i], (∇X1, ∇X2)) = withgradient(X1, X2) do X1, X2
			W1 = sum(X1 ∩ B1) ⊘ sum(X1 ∩ C)
			W2 = sum(X2 ∩ B2) ⊘ sum(X2 ∩ C)
			val(W1 ⊗ W2)
		end

		# The gradient step has to be performed in the real-semiring 
		# (not the log-semiring). `val.(X)` casts `X` in the 
		# real-semiring while `K.(...)` converts it back 
		# to the log-semiring.
		X1 = K.( val.(X1) + (τ ⊗ ∇X1) )
		X2 = K.( val.(X2) + (τ ⊗ ∇X2) )
	end

	plot(F, ylabel="F", xlabel = "step", legend = false)
end

# ╔═╡ 0bb4f2f7-712b-418f-8066-5d523775252f
sum(A2 ∩ B1) ⊘ sum(A2 ∩ C)

# ╔═╡ 9764dfdf-9427-45c0-8381-a16ca6c1e109
X1

# ╔═╡ b3b081ee-87b0-41ca-94f6-4cc4bed8a07d
X2

# ╔═╡ Cell order:
# ╟─08803e3f-cf88-4d91-962f-220709e51019
# ╠═00f12bce-d9dc-11ed-2b99-edb4fd8c8b94
# ╟─1e8e64d2-38be-4474-af8e-077426bbdc92
# ╠═86528707-e1b5-4cbc-bc06-270ead357ada
# ╟─2bde23e4-010f-4102-97b0-5500e1e85fb0
# ╟─ccb356bb-9592-4f04-a4fe-79d7d2e439f4
# ╟─29fb4884-84b5-4e54-a339-50703c1bb89f
# ╟─70b8e8f7-e6c1-4dbc-b165-6eb301395852
# ╟─f7619ff9-7142-4964-b981-d6603bb3df07
# ╟─a09dd6a4-560c-403d-8278-819f16c8451a
# ╟─1c3a8f7e-71e4-4db7-8cc1-9a327d006780
# ╠═2ef2e1c1-6872-4d63-8d15-8b118cbad37f
# ╠═0bb4f2f7-712b-418f-8066-5d523775252f
# ╠═9764dfdf-9427-45c0-8381-a16ca6c1e109
# ╠═b3b081ee-87b0-41ca-94f6-4cc4bed8a07d
