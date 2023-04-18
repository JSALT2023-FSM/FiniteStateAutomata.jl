### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 00f12bce-d9dc-11ed-2b99-edb4fd8c8b94
begin
	using Pkg
	Pkg.activate("../")

	using Revise
	using FiniteStateAutomata
	using ChainRulesCore
	using SparseArrays
	using Semirings
	using Zygote

	using PlutoUI
end

# ╔═╡ 86528707-e1b5-4cbc-bc06-270ead357ada
K = LogSemiring{Float32}

# ╔═╡ 72b2d61b-07ac-4be6-b4f1-e8df9e1ffb12
K(log(2)) / K(log(3))

# ╔═╡ 9d98c862-bb25-4b89-93e8-8d571fe49bcd
K(log(2/3))

# ╔═╡ 099cc966-d043-46f8-8439-cb4f1ca00a6d
function ChainRulesCore.rrule(::typeof(val), x::LogSemiring) 
	function pullback(Δy) 
		(NoTangent(), Δy)
	end
	val(x), pullback
end

# ╔═╡ 573ce00f-30c7-4eac-8ed2-38d1de31907d
function ChainRulesCore.rrule(::Type{<:LogSemiring}, x) 
	function pullback(Δy)
		(NoTangent(), Δy)
	end
	K(x), pullback
end

# ╔═╡ e7993dcf-8a88-4f83-a76c-774af3628bf6
typeof(+)

# ╔═╡ 8a7fe3aa-ba81-4f08-b2b2-c309225ce86d
function ChainRulesCore.rrule(::typeof(*), x::LogSemiring, y::LogSemiring) 
	function pullback(Δy) 
		println("* ", typeof(Δy))
		(NoTangent(), 1, 1)
	end
	x * y, pullback
end

# ╔═╡ bbd5a9b9-9275-4100-be99-cec0a0a1fb39
function add(x::LogSemiring, y::LogSemiring)
	x + y
end

# ╔═╡ 54bff4e6-73ee-4b57-872c-7e6b610b3e0e
function ChainRulesCore.rrule(::typeof(add), x::LogSemiring, y::LogSemiring) 
	z = x + y
	function pullback(Δy) 
		println("+ ", typeof(Δy))
		exp_z = exp(val(z))
		(NoTangent(), Δy * exp(val(x)) / exp_z, Δy * exp(val(y)) / exp_z)
	end
	z, pullback
end

# ╔═╡ e426112e-0a95-4f86-8435-124566becf69
begin	
	gradient(2) do x
		val(add(K(x), K(3)))
	end 
end

# ╔═╡ fbafa614-9124-49cb-9aea-6b34817675eb
typeof(LogSemiring{Float32}) 

# ╔═╡ cccd5ff6-064f-4a51-9076-c1fd4f5aaf2f
Float32.(x)

# ╔═╡ 29fb4884-84b5-4e54-a339-50703c1bb89f
A = FST(
	sparsevec([1, 2], one(K), 4),
	sparse([1, 1, 1, 2, 4], [2, 3, 4, 3, 2], one(K), 4, 4),
	sparsevec([3, 4], one(K), 4),
	zero(K),
	["a", "b", "d", "c"]
)

# ╔═╡ 4e07e8a8-16a6-48e9-b7f6-f1c659672cac
A 

# ╔═╡ 7b0d7d12-0c05-4877-ac02-05050d42acb3
A

# ╔═╡ b68cd52e-b0c5-4bc6-9f15-5e9c8e67a608
A

# ╔═╡ ceaf7631-ca59-4c8a-a4fa-0223da2d9fe7
A

# ╔═╡ ab2853b8-61c2-43f1-9010-9ffbd26457a5
gA = gradient(A) do A
	# C = A ∩ B
	# num = sum(C)
	# D = A ∩ E
	# den = sum(D)
	# loss = num / den
	X = filter(q -> λ(A)[q] != "d", A)
	sum(X)
end[1]

# ╔═╡ 8135a193-5fbc-416f-8475-600f2caf2553


# ╔═╡ ea87c102-7156-41c9-b6d9-61920275eeee
gA

# ╔═╡ 58003de3-f64a-4c90-8ef7-cb5f7c5dd6e4
X = filter(q -> λ(A)[q] != "d", A)

# ╔═╡ 4940ad16-3578-44a5-a241-0f89eedfac55
gradient(A) do A
	# C = A ∩ B
	# num = sum(C)
	# D = A ∩ E
	# den = sum(D)
	# loss = num / den
	X = A
	sum(X)
end[1]

# ╔═╡ ae1183b6-691e-48e7-8fac-350e18fcf806
A 

# ╔═╡ 5073f367-3c37-4f91-ba9c-13520d2a9f98
statemap(A, ["a", "b", "c"]) do q
	λ(A)[q] == "a" && return [1]
	λ(A)[q] == "b" && return [2]
	λ(A)[q] == "c" && return [3]
	return []
end

# ╔═╡ 8e6024ec-1fa1-4e7a-8b37-a44193f048aa
Base.filter(q -> λ(A)[q] != "d", A)

# ╔═╡ c884af39-3785-4c6b-956f-48098bf37dce
sum(A)

# ╔═╡ f7619ff9-7142-4964-b981-d6603bb3df07
B = FST(
	sparsevec([1], one(K), 3),
	sparse([1, 1, 2, 2], [2, 3, 1, 3], one(K), 3, 3),
	sparsevec([2, 3], one(K), 3),
	zero(K),
	["a", "b", "c"]
)

# ╔═╡ 202c4c52-c0b8-43f8-a8bd-a0e26f887c4d
(v1, g1) = withgradient(A) do A
	X = A ∩ B
	Y = A ∩ B
	#one(K) / sum(Y)
	sum(X) / one(K)
end

# ╔═╡ 4c81bbeb-e889-4bc9-af7f-27a18f8a6109
g1[1]

# ╔═╡ b6d14311-8359-403e-b669-e9739b5d913f
withgradient(A ∩ B) do X
	sum(X)
end.grad[1]

# ╔═╡ c7a86d45-509d-4a58-bbac-c3844178edf9
B

# ╔═╡ 1c3a8f7e-71e4-4db7-8cc1-9a327d006780
C = FST(
	sparsevec([1, 3], one(K), 4),
	sparse([1, 1, 2, 2, 3, 4], [2, 3, 3, 4, 4, 1], one(K), 4, 4),
	sparsevec([3, 4], one(K), 4),
	zero(K),
	["a", "b", "c", "d"]
)

# ╔═╡ 5d75ae56-227a-46b3-b92f-5d4792093ee7
(v2, g2) = withgradient(A + g1[1]) do A
	X = A ∩ B
	Y = A ∩ C
	sum(Y) / sum(X) 
end

# ╔═╡ deb0302d-daf9-412e-ab72-05013bba937d
(v3, g3) = withgradient(A + g1[1] + g2[1]) do A
	X = A ∩ B
	Y = A ∩ C
	sum(Y) / sum(X) 
end

# ╔═╡ 5d0c3ed2-8fc1-439f-ab83-07047e6fd931
(v4, g4) = withgradient(A + g1[1] + g2[1] + g3[1]) do A
	X = A ∩ B
	Y = A ∩ C
	sum(Y) / sum(X) 
end

# ╔═╡ 5b454873-e878-4516-ae69-ca857227c0bd
(v5, g5) = withgradient(A + g1[1] + g2[1] + g3[1] + g4[1]) do A
	X = A ∩ B
	Y = A ∩ C
	sum(Y) / sum(X) 
end

# ╔═╡ dec671cc-e6dd-49b5-8c4c-fc11a8eb0450
(v6, g6) = withgradient(A + g5[1]) do A
	X = A ∩ B
	Y = A ∩ C
	sum(Y) / sum(X) 
end

# ╔═╡ a549d176-00b9-411d-9ec2-7fa3762d60c6
A ∩ C

# ╔═╡ d8eed0bd-b75a-4973-81ee-aa739492fdd3
A 

# ╔═╡ c99e23c0-eb5f-4d63-b198-70ab7d73d542
A ∩ C

# ╔═╡ f56f508a-124c-40cf-81c5-e0a06203a7d1
C

# ╔═╡ 56fe0b5a-2b9d-4596-bfe4-87b6e36276a3
A

# ╔═╡ 5dbcaf39-9654-4e10-a67d-2c4cfc6a7240
filter(A) do q
	λ(A)[q] != "d"
end

# ╔═╡ 39f8b673-5154-4af1-9e84-9aacead65095
maximum([3, 2, 1])

# ╔═╡ 9068d3f4-cc9e-412c-b1e3-d9486d45e8d9
filter(q -> λ(A)[q] != "d", A)

# ╔═╡ d25c2762-7c89-4647-8550-8fa62485295d
gradient(filter(q -> λ(A)[q] != "d", A)) do A
	sum(A)
end[1]

# ╔═╡ 3dce0aed-83c3-443f-b3c3-dd7706ea50e4
@bind n Slider(0:(nstates(A)+1))

# ╔═╡ 506ae657-3caa-4bea-9eb0-0d587d6e1ab3


# ╔═╡ 4e923688-9fcf-4d45-9e19-e8c6d7be5539
log(2)

# ╔═╡ 8d4bb635-e1e5-4272-ad09-b2fa2feaace6
val(K(1) + K(2))  |> exp

# ╔═╡ 23dac3f9-f9e3-495f-b506-b779d23e807c
exp(1) + exp(2)

# ╔═╡ a529fbe3-758a-4446-b374-68f3f8e1bead
function nsteps(A::AbstractFST, n)
	n < 1 && return sparsevec([], one(K), nstates(A))
	
	next = iterate(A)
	iter = 1
	while ! isnothing(next) && iter < n
		next = iterate(A, next[2])
		iter += 1
	end
	isnothing(next) ? sparsevec([], one(K), nstates(A)) : next[1]
end

# ╔═╡ 46cc161f-c459-4f7e-92c8-e99636a923a4
activestates(tokens) = findnz(tokens)[1]

# ╔═╡ aeb32bea-570d-4608-9c0b-c68f434fd942
A, nsteps(A , n) |> activestates

# ╔═╡ 975c35e0-d9e8-49c4-975a-19e028e3a8bb
A, nsteps(A |> reverse, n) |> activestates

# ╔═╡ 07435715-fc51-414a-8143-e76da953d3dc
FiniteStateAutomata.dot_write(stdout, A)

# ╔═╡ 3273ec9b-9166-45b2-aad0-70e0b82349f0
A |> reverse 

# ╔═╡ 849f0312-421c-480f-9f1b-864c418a2227
FST(A |> reverse)

# ╔═╡ Cell order:
# ╠═00f12bce-d9dc-11ed-2b99-edb4fd8c8b94
# ╠═86528707-e1b5-4cbc-bc06-270ead357ada
# ╠═72b2d61b-07ac-4be6-b4f1-e8df9e1ffb12
# ╠═9d98c862-bb25-4b89-93e8-8d571fe49bcd
# ╠═099cc966-d043-46f8-8439-cb4f1ca00a6d
# ╠═573ce00f-30c7-4eac-8ed2-38d1de31907d
# ╠═e7993dcf-8a88-4f83-a76c-774af3628bf6
# ╠═8a7fe3aa-ba81-4f08-b2b2-c309225ce86d
# ╠═bbd5a9b9-9275-4100-be99-cec0a0a1fb39
# ╠═54bff4e6-73ee-4b57-872c-7e6b610b3e0e
# ╠═e426112e-0a95-4f86-8435-124566becf69
# ╠═fbafa614-9124-49cb-9aea-6b34817675eb
# ╠═202c4c52-c0b8-43f8-a8bd-a0e26f887c4d
# ╠═4c81bbeb-e889-4bc9-af7f-27a18f8a6109
# ╠═5d75ae56-227a-46b3-b92f-5d4792093ee7
# ╠═deb0302d-daf9-412e-ab72-05013bba937d
# ╠═5d0c3ed2-8fc1-439f-ab83-07047e6fd931
# ╠═5b454873-e878-4516-ae69-ca857227c0bd
# ╠═dec671cc-e6dd-49b5-8c4c-fc11a8eb0450
# ╠═b6d14311-8359-403e-b669-e9739b5d913f
# ╠═a549d176-00b9-411d-9ec2-7fa3762d60c6
# ╠═4e07e8a8-16a6-48e9-b7f6-f1c659672cac
# ╠═7b0d7d12-0c05-4877-ac02-05050d42acb3
# ╠═cccd5ff6-064f-4a51-9076-c1fd4f5aaf2f
# ╠═b68cd52e-b0c5-4bc6-9f15-5e9c8e67a608
# ╠═c7a86d45-509d-4a58-bbac-c3844178edf9
# ╠═29fb4884-84b5-4e54-a339-50703c1bb89f
# ╠═ceaf7631-ca59-4c8a-a4fa-0223da2d9fe7
# ╠═ab2853b8-61c2-43f1-9010-9ffbd26457a5
# ╠═8135a193-5fbc-416f-8475-600f2caf2553
# ╠═ea87c102-7156-41c9-b6d9-61920275eeee
# ╠═58003de3-f64a-4c90-8ef7-cb5f7c5dd6e4
# ╠═4940ad16-3578-44a5-a241-0f89eedfac55
# ╠═ae1183b6-691e-48e7-8fac-350e18fcf806
# ╠═5073f367-3c37-4f91-ba9c-13520d2a9f98
# ╠═8e6024ec-1fa1-4e7a-8b37-a44193f048aa
# ╠═c884af39-3785-4c6b-956f-48098bf37dce
# ╠═f7619ff9-7142-4964-b981-d6603bb3df07
# ╠═1c3a8f7e-71e4-4db7-8cc1-9a327d006780
# ╠═d8eed0bd-b75a-4973-81ee-aa739492fdd3
# ╠═c99e23c0-eb5f-4d63-b198-70ab7d73d542
# ╠═f56f508a-124c-40cf-81c5-e0a06203a7d1
# ╠═56fe0b5a-2b9d-4596-bfe4-87b6e36276a3
# ╠═5dbcaf39-9654-4e10-a67d-2c4cfc6a7240
# ╠═39f8b673-5154-4af1-9e84-9aacead65095
# ╠═9068d3f4-cc9e-412c-b1e3-d9486d45e8d9
# ╠═d25c2762-7c89-4647-8550-8fa62485295d
# ╠═3dce0aed-83c3-443f-b3c3-dd7706ea50e4
# ╠═aeb32bea-570d-4608-9c0b-c68f434fd942
# ╠═975c35e0-d9e8-49c4-975a-19e028e3a8bb
# ╠═506ae657-3caa-4bea-9eb0-0d587d6e1ab3
# ╠═4e923688-9fcf-4d45-9e19-e8c6d7be5539
# ╠═8d4bb635-e1e5-4272-ad09-b2fa2feaace6
# ╠═23dac3f9-f9e3-495f-b506-b779d23e807c
# ╠═a529fbe3-758a-4446-b374-68f3f8e1bead
# ╠═46cc161f-c459-4f7e-92c8-e99636a923a4
# ╠═07435715-fc51-414a-8143-e76da953d3dc
# ╠═3273ec9b-9166-45b2-aad0-70e0b82349f0
# ╠═849f0312-421c-480f-9f1b-864c418a2227
