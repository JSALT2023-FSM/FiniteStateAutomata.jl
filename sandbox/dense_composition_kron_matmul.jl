### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 0a906f46-1ee0-4814-860a-56e127b71593
begin
	using Pkg
	Pkg.develop(path="../../finitestateautomata.jl/")
	
	using Revise
	using FiniteStateAutomata
    using PlutoUI, BenchmarkTools
	using Profile
end

# ╔═╡ 8b30c1e2-f00d-4d8c-b91a-9d7401a4628f
S = TropicalSemiring{Float32}

# ╔═╡ 7e5488a5-4de3-4543-91f5-09e8db4f5c36
symbols = Dict(1 => "a", 2 => "b", 3 => "c")

# ╔═╡ 08e20858-4f43-4614-9d64-fcf7a32c7fa8
function dense_composition_matmul(A,B)
	S = semiring(A)
	MA, MB = M(A), M(B)

	n = numstates(A)
	m = numstates(B)

	ml = max(size(MA, 4), size(MB, 3))
	
	MA = cat(MA, 
		zeros(S, n, n, size(MA, 3), ml - size(MA, 4)); dims=4)
	MB = cat(MB, 
		zeros(S, m, m, ml - size(MB, 3), size(MB, 4) ); dims=3)
	MC = zeros(S,(n*m,n*m, size(MA,3), size(MB,4)))
	for i in 1:n
		for j in 1:n
			for k in 1:m
				for l in 1:m
					MC[(k-1)*n+i,(l-1)*n+j,:,:] .= MA[i,j,:,:]*MB[k,l,:,:]
				end
			end
		end
	end
	TensorFST(MC, kron(α(A), α(B)), kron(ω(A), ω(B)))
end

# ╔═╡ b3549577-dba3-4907-9e65-917918b02a97
function dense_composition_kron(fstA, fstB)	
	S = semiring(fstA)
	MA, MB = M(fstA), M(fstB)

	pMA = permutedims(MA, (3, 4, 1, 2))
	pMB = permutedims(MB, (3, 4, 1, 2))

	ml = max(size(pMA, 2), size(pMB, 1))
	pMA = cat(pMA, 
		zeros(S, size(pMA, 1), ml - size(pMA, 2),numstates(fstA), numstates(fstA)); dims=2)
	pMB = cat(pMB, 
		zeros(S, ml - size(pMB, 1), size(pMB, 2), numstates(fstB), numstates(fstB)); dims=1)
	
	Q = numstates(fstA) * numstates(fstB)
	MC = zeros(S, size(pMA, 1), size(pMB, 2), Q, Q)
	for x in 1:size(MC, 1)
		for z in 1:size(MC, 2)
				MC[x,z,:,:] = sum(
					y -> kron(pMA[x,y,:,:], pMB[y,z,:,:]), 
					1:size(pMA, 2) ; 
					init=zeros(S, Q, Q)
				)
		end
	end

	TensorFST(permutedims(MC, (3, 4, 1, 2)), kron(α(fstA), α(fstB)), kron(ω(fstA), ω(fstB)))
end

# ╔═╡ 8cb74c7d-e329-4b26-9286-4014a1db51dc
A = VectorFST(
	[
		Arc{S}[(2, 1, 2, S(0.1)), (3, 2, 1, S(0.2))],
		Arc{S}[(2, 3, 1, S(0.3)), (4, 1, 1, S(0.4))],
		Arc{S}[(4, 2, 2, S(0.5))],
		Arc{S}[]
	],
	1,
	[zero(S),zero(S),zero(S),S(0.6)]
)

# ╔═╡ dff11a7c-6bee-4611-842b-a8f070c2aa0a
draw(A; symbols)

# ╔═╡ ce4f48c1-e6ca-4ce1-948c-9192020c9f95
B = VectorFST(
	[
		Arc{S}[(2, 2, 3, S(0.3))],
		Arc{S}[(3, 1, 2, S(0.4))],
		Arc{S}[(3, 1, 2, S(0.6))],
	],
	1,
	S[zero(S),zero(S),S(0.7)]
)

# ╔═╡ 9e46b70d-f904-43b5-9f1a-3a6f52402c2d
draw(B; symbols)

# ╔═╡ ce623bef-6f01-483c-b262-54d60e7751dc
begin
	fstA = convert(TensorFST{S, Array{S,4}}, A)
	draw(fstA, symbols=symbols)
end

# ╔═╡ f0af788f-adb9-440c-b6b2-d0c5d8ad0943
begin
	fstB = convert(TensorFST{S, Array{S,4}}, B)
	draw(fstB, symbols=symbols)
end

# ╔═╡ 679bb9f8-801a-4380-ab20-02af120f246e
md"Now with TensorFST"

# ╔═╡ 01b4e07a-ab18-45c3-b3ba-3ad4be2a101d
fstC = fst_composition(fstA, fstB, length(symbols))

# ╔═╡ 9972ca24-787a-48dc-9032-de3516b3a368
draw(fstC; symbols)

# ╔═╡ fc385050-8951-40be-a42f-5e870ebd0d4b


# ╔═╡ 835b438d-5496-42e3-9613-43ddad812802
C = dense_composition_kron(fstA, fstB)

# ╔═╡ 795ffba4-537c-4a31-90af-1f53794e075b
draw(C; symbols)

# ╔═╡ 4562b52c-9c6d-475e-a569-bbf0c68420e0
C2 = dense_composition_matmul(fstA, fstB)

# ╔═╡ d1143422-d9d6-4102-aa00-0e317a447687
draw(C2; symbols)

# ╔═╡ 357c986e-3b01-4751-8b86-4aa00b04f8c3
@benchmark dense_composition_matmul(fstA, fstB)

# ╔═╡ b13d4338-6987-444f-a023-dc43a55a35fb
@benchmark dense_composition_kron(fstA, fstB)

# ╔═╡ Cell order:
# ╠═0a906f46-1ee0-4814-860a-56e127b71593
# ╠═8b30c1e2-f00d-4d8c-b91a-9d7401a4628f
# ╠═7e5488a5-4de3-4543-91f5-09e8db4f5c36
# ╠═08e20858-4f43-4614-9d64-fcf7a32c7fa8
# ╠═b3549577-dba3-4907-9e65-917918b02a97
# ╠═8cb74c7d-e329-4b26-9286-4014a1db51dc
# ╠═dff11a7c-6bee-4611-842b-a8f070c2aa0a
# ╠═ce4f48c1-e6ca-4ce1-948c-9192020c9f95
# ╠═9e46b70d-f904-43b5-9f1a-3a6f52402c2d
# ╠═ce623bef-6f01-483c-b262-54d60e7751dc
# ╠═f0af788f-adb9-440c-b6b2-d0c5d8ad0943
# ╟─679bb9f8-801a-4380-ab20-02af120f246e
# ╠═01b4e07a-ab18-45c3-b3ba-3ad4be2a101d
# ╠═9972ca24-787a-48dc-9032-de3516b3a368
# ╠═fc385050-8951-40be-a42f-5e870ebd0d4b
# ╠═835b438d-5496-42e3-9613-43ddad812802
# ╠═795ffba4-537c-4a31-90af-1f53794e075b
# ╠═4562b52c-9c6d-475e-a569-bbf0c68420e0
# ╠═d1143422-d9d6-4102-aa00-0e317a447687
# ╠═357c986e-3b01-4751-8b86-4aa00b04f8c3
# ╠═b13d4338-6987-444f-a023-dc43a55a35fb
