### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ b4e4c248-da19-11ed-160a-a708435cbbff
begin
	import Pkg
	Pkg.activate("./mmi")
	using Revise
	using Semirings
	using SparseArrays
	using LinearAlgebra
	using FiniteStateAutomata
end

# ╔═╡ 9418be4b-0639-4715-bfde-ae8f3f119a8c
begin
	K = LogSemiring{Float32}
	Σ = ["a", "b"]
	λ = repeat(Σ, 3)
	N = 3
	H = (K∘log).(rand(length(Σ), N+1))
end

# ╔═╡ 8c5250b6-51b8-4b89-b506-38a4ecf93f86
A = DenseFSA(
	# H,
	K.(H),
	Σ,
	zero(K),
)

# ╔═╡ ca6b6f54-1082-4081-9bbf-0ec46106a0cf
B = FSA(
	sparsevec([1, 2], one(K), length(λ)),
	sparse([1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6], [3, 4, 5, 6, 3, 4, 5, 6, 3, 4, 5, 6], one(K), length(λ), length(λ)),
	sparsevec(1:length(λ), one(K), length(λ)),
	one(K),
	λ,
)

# ╔═╡ 30492b65-4a36-44ad-a68c-19d022f6cf11
C = A ∩ B

# ╔═╡ 131ee7b4-9b7a-4833-989f-b367070f7544
md"""
# Test `sum`
"""

# ╔═╡ db55b5fa-5829-4612-98d2-d95432fdd31e
begin
	G = C.C' * H
	tmp = T(B) * (G[:, 2] .* (T(B) * (G[:, 3] .* G[:, 4] .* ω(B))))
	W = dot(α(B) .* G[:, 1], tmp)
end

# ╔═╡ 0f2a4a60-fe4a-416d-ac68-19cac33c7d3c
sum(C)

# ╔═╡ c797063d-a347-4f0d-bc73-1617337948d3
md"### Forward"

# ╔═╡ 70ed6238-6b0c-494d-b1f7-acbd4eadbcf6
begin
	u = fill!(similar(G, size(G, 1), N), zero(K))
	u[:, 1] = α(B)
	for n in 2:N
		u[:, n] = T(B)' * (u[:, n - 1] .* G[:, n - 1])
	end
	u
end

# ╔═╡ d6d59daf-7c2e-46d1-87fa-d77b9eeba67f
md"""
### Backward
"""

# ╔═╡ 77d58810-ff50-406d-8a9c-d09ffbabd78a
begin
	v = fill!(similar(G, size(G, 1), N), zero(K))
	v[:, end] = G[:, end] .* ω(B)
	for n in N - 1: -1: 1
		v[:, n] = T(B) * (G[:, n+1] .* v[:, n+1])
	end
	v
end

# ╔═╡ a5b93f32-a554-471f-a1f1-b74ec4ccec60
md"### Test"

# ╔═╡ 3424db3c-5f84-440d-9b7d-809d90ecd247
sum(C),
dot(G[:, end - 1] .* u[:, end], G[:, end] .* ω(B)),
dot(G[:, 1] .* α(B), v[:, 1])

# ╔═╡ 468647c0-9baf-40c9-9599-442057a6b5f4
begin
	W2 = similar(G, N)
	for n in 1:N
		W2[n] = dot(u[:, n] .* G[:, n], v[:, n])
	end

	W3 = similar(G, N)
	for n in 1:N
		W3[n] = dot(u[:, n], v[:, n] .* G[:, n])
	end
	sum(C), W2, W3
end

# ╔═╡ a593f88f-1b65-4f2c-8d5a-eafff1f2b174
md"## Gradient"

# ╔═╡ 2da4592a-8b2e-41b9-b2b4-42a1bcb28b3a
begin
	∇ₕW = fill!(similar(H), zero(eltype(H)))
	∇ₕW[:, 1:end-1] = C.C * (u .* v)
	∇ₕW[:, end] = C.C * ((C.C' * C.A.H[:, end - 1]) .* u[:, end] .* ω(C.B))
	∇ₕW
end

# ╔═╡ f854585b-c16a-4cbc-b161-dfc44beb583f
gradient(sum, C)

# ╔═╡ Cell order:
# ╠═b4e4c248-da19-11ed-160a-a708435cbbff
# ╠═9418be4b-0639-4715-bfde-ae8f3f119a8c
# ╠═8c5250b6-51b8-4b89-b506-38a4ecf93f86
# ╠═ca6b6f54-1082-4081-9bbf-0ec46106a0cf
# ╠═30492b65-4a36-44ad-a68c-19d022f6cf11
# ╠═131ee7b4-9b7a-4833-989f-b367070f7544
# ╠═db55b5fa-5829-4612-98d2-d95432fdd31e
# ╠═0f2a4a60-fe4a-416d-ac68-19cac33c7d3c
# ╠═c797063d-a347-4f0d-bc73-1617337948d3
# ╠═70ed6238-6b0c-494d-b1f7-acbd4eadbcf6
# ╠═d6d59daf-7c2e-46d1-87fa-d77b9eeba67f
# ╠═77d58810-ff50-406d-8a9c-d09ffbabd78a
# ╠═a5b93f32-a554-471f-a1f1-b74ec4ccec60
# ╠═3424db3c-5f84-440d-9b7d-809d90ecd247
# ╠═468647c0-9baf-40c9-9599-442057a6b5f4
# ╠═a593f88f-1b65-4f2c-8d5a-eafff1f2b174
# ╠═2da4592a-8b2e-41b9-b2b4-42a1bcb28b3a
# ╠═f854585b-c16a-4cbc-b161-dfc44beb583f
