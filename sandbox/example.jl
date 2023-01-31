### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ e783cfcc-8b82-11ed-147f-d18598bc60ff
begin
	using Pkg
	Pkg.activate("../")

	using Revise
	using FSALinalg
	using LinearAlgebra
	using SparseArrays
	using Semirings
end

# ╔═╡ 33a02e4e-44a6-4d66-9a1a-0372b96695d6
md"""
# Using Linear Algebra to manipulate Finite State Automata
"""

# ╔═╡ 0f68bc4e-3fcb-4d28-bba2-4d0032c80981
md""" 
## Matrix-based representation of FSA

We define a FSA over a semiring ``K`` as a 6-tuple: 
```math
\mathcal{M} = (\boldsymbol{\alpha}, \mathbf{T}, \boldsymbol{\omega}, \rho, \boldsymbol{\lambda})
```
where:
* ``\boldsymbol{\alpha} \in K^D`` is the vector of initial states 
* ``\mathbf{T} \in K^{D \times D}`` is the adjacency matrix of the FSA.
* ``\boldsymbol{\omega} \in K^D`` is the vector of final states
* ``\rho \in K`` is the weight associated with the emtpy string ``\epsilon``
* ``\boldsymbol{\lambda} \in K^D`` is the vector of final states

"""

# ╔═╡ 87e1cf8e-da52-4c85-ab6a-bb352a84d779
K = ProbSemiring{Float32}

# ╔═╡ 8602a2b5-16ca-4cec-8041-c6362aa66a48
M1 = FSA(
	sparsevec([3, 2], K[0.5, 1.5], 3),
	sparse([3, 2, 2], [1, 1, 3], K[2.5, 2.5, 1.5], 3, 3),
	sparsevec([1], K[3.5], 3),
	K(0.5),
	["a", "b", "c"]
)

# ╔═╡ 3e34e635-3565-437d-acbe-6bd170bf4382
findall([true, false, true])

# ╔═╡ 028ff41d-f9a7-4ff8-b254-7399be161a48
Diagonal(Symmetric(randn(2, 2)))

# ╔═╡ fe7c8892-2439-4b1d-b661-92f8bb6825bc
StochasticFSA(M1)

# ╔═╡ 9400b791-605f-4cb5-a3d7-71af5a034504
M2 = FSA(
	sparsevec([3, 2], K[0.2, 3], 4),
	sparse([3, 3, 2, 2, 1], [2, 1, 2, 4, 4], K[0.1, 0.2, 2, 0.3, 0.4], 4, 4),
	sparsevec([4], K[0.6], 4),
	zero(K),
	["c", "b", "a", "d"]
)

# ╔═╡ 78aa5aa2-87a7-4bf6-80a2-71fd79bf03d1
parent(M1)

# ╔═╡ 45d44224-9e88-482c-95d8-d1a793b3fd79
renorm(M1)

# ╔═╡ 6d6c0896-dbd3-41d1-8f90-04120bd26087
cumsum(M1 |> renorm)

# ╔═╡ 1bfd0d49-8833-4b05-978c-9830d5c6bc33
M1

# ╔═╡ e7a6ab1c-262e-4a3b-a34b-e2485dfbeedf
R = ProductSemiring{Tuple{ProbSemiring{Float64}, UnionConcatSemiring{StringMonoid}}}

# ╔═╡ 992b92d0-e402-4cec-a232-0a2c6e71d9c7
one(zero(R))

# ╔═╡ 57c74af2-73fe-4c57-af46-7eedff77e119
T1 = UnionConcatSemiring{StringMonoid}

# ╔═╡ edf219c8-9663-4af9-8ace-3ace909fad07
T1(Set(StringMonoid["b"])) * T1(Set(StringMonoid["a"]))

# ╔═╡ ff37939d-195b-4a80-863d-89baf53c2341
convert(M1) do v, l
	if isnothing(l) && iszero(v)
		return ProductSemiring(
			(
				v, 
				zero(UnionConcatSemiring{StringMonoid})
			)
		)
	elseif isnothing(l)
		return ProductSemiring(
			(
				v, 
				one(UnionConcatSemiring{StringMonoid})
			)
		)
	else
		return ProductSemiring((v, UnionConcatSemiring{StringMonoid}(Set(StringMonoid[l]))))
	end
end

# ╔═╡ 08ad3fe4-5be2-4512-b11e-07c6c0a9b786
function Base.isequal(A1::AbstractFSA{K}, A2::AbstractFSA{K}) where K
	L = UnionConcatSemiring{StringMonoid}
	R = ProductSemiring{Tuple{K, L}}
	

	Dₗ = spdiagm(
		[R((one(K), Set(StringMonoid[l]))) for l in A1.λ]
	)
	[ProductSemiring((v, L(Set(StringMonoid[λ(A1)[i]])))) for (i, v) in zip(initstates(A1)...)]
	[ProductSemiring((v, L(Set(StringMonoid[λ(A1)[i]])))) for (i, v) in zip(finalstates(A1)...)]
	# FSA(
	# 	Dₗ * sparsevec(initstates(A1)[1], one())
	# )
end

# ╔═╡ 6b4504c0-6775-49bf-91fb-759c8d528395
isequal(M1, M2)

# ╔═╡ 7a4e1c92-e01a-4b08-9b98-416a502aef17
L = [R((one(ProbSemiring{Float64}), Set(StringMonoid[l]))) for l in M1.λ]

# ╔═╡ b34b0f0c-1bfc-4e4b-a388-8550fd62d628
SM1 = FSA(
	spdiagm(L) * sparsevec(findnz(α(M1))[1], one(R), nstates(M1)),
	sparse(findnz(T(M1))[1:2]..., one(R), nstates(M1), nstates(M1)) * spdiagm(L),
	sparsevec(findnz(ω(M1))[1], one(R), nstates(M1)),
	! iszero(ρ(M1)) ? one(R) : zero(R),
	M1.λ
)

# ╔═╡ a7b8be4a-a77c-4b26-85d8-f38b2d91965d
cumsum(SM1) == cumsum(SM1)

# ╔═╡ 780eb548-9bb3-40a4-8c95-f5dbb371f070
sum(T(SM1)' * spdiagm(α(SM1)), dims = 2)

# ╔═╡ 6413cd1f-22d1-4cc0-bd03-39d3df88c823
sum(spdiagm(α(SM1)) * T(SM1), dims=2)

# ╔═╡ 0fae1d09-a501-43d3-89f7-2245f6f6d1fd
T(SM1)' * α(SM1)

# ╔═╡ a7e8863d-73cb-4921-8c1a-cd68ba3c5e5b
v₀ = α(SM1) .* λ(SM1) 

# ╔═╡ ce204f13-d907-47df-84bc-7a6aca678583
SM1.T' * SM1.T' * v₀

# ╔═╡ 50a6809e-c10d-41be-9217-00512fd387e2
md"""
## FSA Operations
"""

# ╔═╡ e252d8ea-9763-4030-9440-59e3a63ec4c2
md"""
### Union

```math
\mathcal{M}_3 = \mathcal{M}_1 \cup \mathcal{M}_2
```
where: 
```math
\begin{align}
\boldsymbol{\alpha}_3 &= \begin{bmatrix}
	\boldsymbol{\alpha}_1 \\
	\boldsymbol{\alpha}_2
\end{bmatrix} & 
\mathbf{T}_3 &= \begin{bmatrix}
	\mathbf{T}_1 & \\
	& \mathbf{T}_2
\end{bmatrix} \\
\boldsymbol{\omega}_3 &= \begin{bmatrix}
	\boldsymbol{\omega}_1 \\
	\boldsymbol{\omega}_2
\end{bmatrix} &
\boldsymbol{\lambda}_3 &= \begin{bmatrix}
	\boldsymbol{\lambda}_1 \\
	\boldsymbol{\lambda}_2
\end{bmatrix} \\
\rho_3 &= \rho_1 \oplus \rho_2 \\
\end{align}
```
"""

# ╔═╡ 590ea7ba-eeae-4f03-8189-a199714001fe
M1 ∪ M2 |> renorm

# ╔═╡ 5a2320c0-69c9-4ab2-b350-767b9db2d228
md"""
### Concatenation

```math
\mathcal{M}_3 = \mathcal{M}_1 \cdot \mathcal{M}_2
```
where: 
```math
\begin{align}
\boldsymbol{\alpha}_3 &= \begin{bmatrix}
	\boldsymbol{\alpha}_1 \\
	\rho_1 \otimes \boldsymbol{\alpha}_2
\end{bmatrix} & 
\mathbf{T}_3 &= \begin{bmatrix}
	\mathbf{T}_1 & \boldsymbol{\omega}_1 \boldsymbol{\alpha}_2^\top \\
	& \mathbf{T}_2
\end{bmatrix} \\
\boldsymbol{\omega}_3 &= \begin{bmatrix}
	\bar{0} \otimes \boldsymbol{\omega}_1 \\
	\boldsymbol{\omega}_2
\end{bmatrix} &
\boldsymbol{\lambda}_3 &= \begin{bmatrix}
	\boldsymbol{\lambda}_1 \\
	\boldsymbol{\lambda}_2
\end{bmatrix} \\
\rho_3 &= \rho_1 \otimes \rho_2 \\
\end{align}
```

"""

# ╔═╡ d89b13de-0673-4cce-ab68-c3ae7bc89d5d
cat(M1, M2)

# ╔═╡ b7a00261-686e-4c29-a11f-42764c60762d
md"""

## Kleene closure

```math
\mathcal{M}_2 = \mathcal{M}_1^*
```
where: 
```math
\begin{align}
\boldsymbol{\alpha}_2 &= \boldsymbol{\alpha}_1 & 
\mathbf{T}_2 &= \mathbf{T}_1 + \boldsymbol{\omega}_1 \boldsymbol{\alpha}_1^\top \\
\boldsymbol{\omega}_2 &= \boldsymbol{\omega}_1 &
\boldsymbol{\lambda}_2 &= \boldsymbol{\lambda}_1 \\
\rho_2 &= \begin{cases} 
	\bar{1} & \text{if } \rho_1 = \bar{0} \\
	\rho_1 & \text{otherwise}
\end{cases}
\end{align}
```
"""

# ╔═╡ 309eb417-5487-4362-990b-6221cdb8fc88
M1

# ╔═╡ 97b63bdc-8655-4713-9908-d04ccc801932
closure(M1)

# ╔═╡ 7299f93d-9318-4778-812f-1cb9246c2110
M2

# ╔═╡ 64083d10-8846-4d18-8893-ff5d04a3941a
closure(M2)

# ╔═╡ 0fcc8354-8139-4c72-acc1-c9378939e275
md"""
## Reversal

```math
\mathcal{M}_2 = \mathcal{M}_1^\top
```
where: 
```math
\begin{align}
\boldsymbol{\alpha}_2 &= \boldsymbol{\omega}_1 & 
\mathbf{T}_2 &= \mathbf{T}_1^\top \\
\boldsymbol{\omega}_2 &= \boldsymbol{\alpha}_1 &
\boldsymbol{\lambda}_2 &= \boldsymbol{\lambda}_1   \\
\rho_2 &= \rho_1 \\
\end{align}
```

"""

# ╔═╡ 24f2c0dc-d99b-4809-9080-103a48eb0083
reverse(M1)

# ╔═╡ cf11045e-c8ef-48f0-8ab7-eb38b5119398
reverse(M2)

# ╔═╡ b63c5e21-2218-42ba-80da-dcd74d0ddc38
md"""
## Weight pushing
"""

# ╔═╡ 4e620b69-903b-4085-b45a-cbd8f12a5788
M1

# ╔═╡ 59b01860-6c09-43d7-ae32-70561fec0e58
sum(M1.T, dims=1)

# ╔═╡ 4977c8c7-8a43-4ef0-ad18-2aaa9c147706
function graph_rowspace(T::AbstractMatrix)
	sum(T, dims=1)
	roots = findall(!iszero, sum(T, dims=1))
	I = getindex.(roots, 2)
	sparse(I, I, one(eltype(T)), size(T, 1), size(T, 2)), I
end

# ╔═╡ 998c5d5a-a28f-43ac-900f-508da56e9d3b
function topologicalorder(M)
	T = M.T
	visited = Set()
	Q = size(T, 1)
	iter = 1
	
	s = sum(M.T, dims = 1)
	roots = getindex.(findall(iszero, s), 2)
	L = roots
	while length(roots) > length(visited) && iter < Q
		# Mark the root vertices as visited.
		visited = visited ∪ Set(roots)
		
		# Compute the standard basis of the rowspace
		I = getindex.(findall(!iszero, s), 2)
		R = sparse(I, I, one(eltype(T)), Q, Q)

		# Project the graph onto the rowspace
		T = R * T
		
		s = sum(T, dims = 1)
		roots = getindex.(findall(iszero, s), 2)
		L = vcat(L, filter(x -> x ∉ visited, roots))

		iter += 1
	end

	nnz(T) > 0 ? throw(ArgumentError("FSA is not acyclic")) : L
end

# ╔═╡ 09a701f0-ea44-419f-9de4-f5b4682e9d94
[1, 2, 3][[2, 1, 3]] 

# ╔═╡ 8a653fd5-f5c5-4fc6-88d7-6eeb4da81216
topologicalorder(M1)

# ╔═╡ 5a3f44fc-2016-49a9-9e32-30c9a2bff081
R * M1.T

# ╔═╡ 67e55bab-f2cd-4181-b9a0-33b2181e207b
roots = findall(iszero, sum(M1.T, dims=1))

# ╔═╡ d3256762-6488-4b11-8490-261a98fad0b2
getindex.(roots, 2)

# ╔═╡ 747e8dbb-d41a-442d-a483-bef308ef7f40
M1.T' * spdiagm(M1.α)

# ╔═╡ e7dd85e6-e03e-48ad-96fd-6b9e7c4cb968
M2

# ╔═╡ 6d64bbd1-dbe5-463d-9fd5-bf72fb92ad45
Tuple{LogSemiring{Float64}, BoolSemiring}.parameters

# ╔═╡ 235dfa07-9c56-47a2-a4a8-0b540625bc26


# ╔═╡ 941c7e7a-6b6d-48a8-8426-426b9ddcb9e1
isbitstype(Core.SimpleVector)

# ╔═╡ 7ad2f7f3-e6b6-4707-ac81-22c51d4a8316
function DFS(M, v, visited, order)
	ws = v == 0 ? M.α : M.T[v,:]
	for w in findnz(ws)[1]
		visited[w] != 0 && continue
		visited[w] = order + 1
		DFS(M, w, visited, order + 1)
	end
	visited
end

# ╔═╡ 6b705783-890d-49f2-b3a2-6228897e2fb4
DFS(M2, 0, [0, 0, 0, 0], 0)

# ╔═╡ 2a1daba5-9035-4a5c-b32d-c73c022135c6
M1

# ╔═╡ 27b79f65-fd20-4970-927c-1fe5a453487c
M1.T

# ╔═╡ dbc9c3c2-e6c6-4d3c-81f0-3b8c3dcd4c02
M1.T'

# ╔═╡ e6ed275d-93b4-4e9b-a05b-5d59915536b6
spdiagm(M1.α)

# ╔═╡ 97fa01d5-ecd4-453d-9963-9b7d50eb67cc
M1.T * M1.T # *M1.T * M1.T

# ╔═╡ 38baae5d-3d15-4c79-840a-59635caf7636
function push(fsa)
	σ = zeros(eltype(fsa.α), length(fsa.α)) 
	visited = zeros(Bool, length(fsa.α))

	vₙ = fsa.α
	σ .+= vₙ
	I, V = findnz(vₙ)
	visited[I] .= true
	println(Array(vₙ))
	while nnz(vₙ) > 0
		vₙ = fsa.T' * vₙ
		σ .+= vₙ
		println(Array(vₙ))
		
		I, V = findnz(vₙ)
		J = findall(i -> ! visited[i], I)
		vₙ = sparsevec(I[J], V[J], length(vₙ))
		
		visited[I] .= true
	end

	σ, visited
end

# ╔═╡ 8915f5e6-e9cd-4ba4-a073-bdb1cb1f7890
push(M2)

# ╔═╡ 5bc565c6-25fd-44e1-8222-3741d0e13afb
length(M1.α)

# ╔═╡ fb468dfe-4176-4b84-90a0-6346fe4c7475
md"""
# Determinization
"""

# ╔═╡ 3b828c13-737a-4b5f-8f1f-02d4c7384b6a
M3 = FSA(
	sparsevec([1, 2, 4], K[1, 1, 1], 6),
	sparse([1, 2, 4, 5, 5, 6, 6], [2, 3, 5, 1, 6, 1, 3], K[1, 1, 1, 1, 1, 1, 1], 6, 6),
	sparsevec([3, 6], K[1, 1], 6),
	one(K),
	["a", "b", "a", "a", "b", "c"]
) 

# ╔═╡ 3b04bdf7-24d8-4269-8576-02ac2bad6c97
C = sparse(
	[1, 2, 3, 4, 5, 6],
	[1, 2, 1, 1, 2, 3],
	1,
	6, 3
)

# ╔═╡ 88ef84f4-b1a6-4159-bbbf-f4f506859950
sort(collect(Set([1, 2, 3])))

# ╔═╡ b2b15ac7-1574-4999-bbc9-3536d3aa18b1
function spmap(K, mapping)
	unique_ids = Set(mapping)
	id2idx = Dict(id => i for (i, id) in enumerate(sort(collect(unique_ids))))
	idx_mapping = [id2idx[id] for id in mapping]
	
	Q, D = length(mapping), length(unique_ids)
	sparse(1:Q, idx_mapping, one(K), Q, D)
end

# ╔═╡ 09ca7aa2-70a7-4d1b-ae9f-bdd937bf5c94
spmap(K, M3.λ)

# ╔═╡ e0e8eb60-7ad7-4015-b226-56af6d405fad
function determinize(M::FSALinalg.AbstractFSA{K}) where K
	C = spmap(K, M.λ)

	# Edges of the determinized FSA.
	edges = Dict()

	# Visited states of the determinized FSA.
	visited = Set()

	# Initialize the stack, `()` represents the initial state 
	# of the determinized FSA.
	stack = Any[((), M.α)]
	while ! isempty(stack)
		# Get the first element from the stack. 
		# `a` is the ancestor state (of the new FSA)
		# and `v` is the children of `a` in `M`. 
		a, v = pop!(stack)

		# Merge states with the same label.
		V = spdiagm(v) * C

		# Iterate over the column of V that have non-zero elements.
		for i in 1:size(V, 2)
			c = V[:, i]
			nnz(c) == 0 && continue

			# A new state is defined based on the indices
			# of the non-zero elements.
			q = Tuple(findnz(c)[1])

			# Record the edge between the ancestor and the new 
			# state.
			edges[a] = push!(get(edges, a, []), q)

			# If the new state `q` has not been visited yet
			# add its children (represented as a vector)
			# on the stack.
			if q ∉ visited 
				push!(stack, (q, M3.T' * c))
				push!(visited, q)
			end
		end
	end
	edges, visited
end

# ╔═╡ 143d6069-6978-43ed-815f-5cf19de83fbf
function fsadet(M, edges, states)
	state2idx = Dict(q => i for (i, q) in enumerate(sort(collect(states))))
	D, Q = length(states), length(M.λ)
	K = eltype(M.α)

	I, V = [], K[]
	for i in edges[()]
		eᵢ = sparsevec(collect(i), one(K), Q)
		push!(I, state2idx[i])
		push!(V, dot(M.α, eᵢ))
	end
	α = sparsevec(I, V, D)
	
	I, J, V = [], [], K[]
	for i in keys(edges)
		i == () && continue
		for j in edges[i]
			push!(I, state2idx[i])
			push!(J, state2idx[j])
			eᵢ = sparsevec(collect(i), one(K), Q)
			eⱼ = sparsevec(collect(j), one(K), Q)
			push!(V, eᵢ' * M.T * eⱼ)
		end
	end
	T = sparse(I, J, V, D, D)

	I, V = [], K[]
	for i in states
		eᵢ = sparsevec(collect(i), one(K), Q)
		push!(I, state2idx[i])
		push!(V, dot(M.ω, eᵢ))
	end
	ω = sparsevec(I, V, D)

	λ = [M.λ[q[1]] for q in sort(collect(states))]

	FSA(α, T, ω, M.ρ, λ)
end

# ╔═╡ 5433be8c-61fa-4262-839a-d9a137c05377
edges, states = determinize(M3)

# ╔═╡ 03289f26-13cf-418f-853a-228b6427cf33
function conv(f::Function, A::AbstractFSA{K}) where K

	sparsevec(initstates(A)[1], [f(i, v) for (i, v) in zip(initstates(A)...)], nstates(A))
	sparse(edges(A)[1], edges(A)[2], [f(j, v) for (i, j, v) in zip(edges(A)...)], nstates(A), nstates(A))
end

# ╔═╡ aea78dcc-02a2-44d4-a24e-acb088d61bc3
edges

# ╔═╡ c8d7966a-61c1-4ab2-856c-6ca5d10648cd
fsadet(M3 |> renorm, edges, states) |> renorm

# ╔═╡ 358d297d-4def-44dc-b79d-6d23aa67d179
sum(M3 |> renorm)

# ╔═╡ c31ccfdc-3687-46cb-9eed-03791bde56d4
sparsevec([2, 3], 1, 5)

# ╔═╡ 75c4a277-9edb-4f47-bd5e-13dcb498bd96
M3 

# ╔═╡ 68b7380b-2514-4362-9504-552fb4e467f7
Z1 = spdiagm(M3.α) * C

# ╔═╡ f2a127b3-df43-486a-a34b-e289539a22bf
for i in 1:size(C, 2)
	v = Z1[:, i]
	println(Tuple(findnz(v)[1]))
end

# ╔═╡ f5806376-9ee8-4991-a84e-d650a149edd1
Z2 = M3.T' * Z1

# ╔═╡ 7f0c719d-fbdf-455a-94f4-a1e899462c5d
Z3 = M3.T' * Z2

# ╔═╡ a8733344-6aba-47a5-bbc4-444fac026d40
C .* Z3

# ╔═╡ Cell order:
# ╟─33a02e4e-44a6-4d66-9a1a-0372b96695d6
# ╠═e783cfcc-8b82-11ed-147f-d18598bc60ff
# ╟─0f68bc4e-3fcb-4d28-bba2-4d0032c80981
# ╠═87e1cf8e-da52-4c85-ab6a-bb352a84d779
# ╠═8602a2b5-16ca-4cec-8041-c6362aa66a48
# ╠═3e34e635-3565-437d-acbe-6bd170bf4382
# ╠═028ff41d-f9a7-4ff8-b254-7399be161a48
# ╠═fe7c8892-2439-4b1d-b661-92f8bb6825bc
# ╠═9400b791-605f-4cb5-a3d7-71af5a034504
# ╠═78aa5aa2-87a7-4bf6-80a2-71fd79bf03d1
# ╠═45d44224-9e88-482c-95d8-d1a793b3fd79
# ╠═6d6c0896-dbd3-41d1-8f90-04120bd26087
# ╠═1bfd0d49-8833-4b05-978c-9830d5c6bc33
# ╠═e7a6ab1c-262e-4a3b-a34b-e2485dfbeedf
# ╠═992b92d0-e402-4cec-a232-0a2c6e71d9c7
# ╠═57c74af2-73fe-4c57-af46-7eedff77e119
# ╠═edf219c8-9663-4af9-8ace-3ace909fad07
# ╠═b34b0f0c-1bfc-4e4b-a388-8550fd62d628
# ╠═03289f26-13cf-418f-853a-228b6427cf33
# ╠═aea78dcc-02a2-44d4-a24e-acb088d61bc3
# ╠═ff37939d-195b-4a80-863d-89baf53c2341
# ╠═08ad3fe4-5be2-4512-b11e-07c6c0a9b786
# ╠═6b4504c0-6775-49bf-91fb-759c8d528395
# ╠═7a4e1c92-e01a-4b08-9b98-416a502aef17
# ╠═a7b8be4a-a77c-4b26-85d8-f38b2d91965d
# ╠═780eb548-9bb3-40a4-8c95-f5dbb371f070
# ╠═6413cd1f-22d1-4cc0-bd03-39d3df88c823
# ╠═0fae1d09-a501-43d3-89f7-2245f6f6d1fd
# ╠═a7e8863d-73cb-4921-8c1a-cd68ba3c5e5b
# ╠═ce204f13-d907-47df-84bc-7a6aca678583
# ╟─50a6809e-c10d-41be-9217-00512fd387e2
# ╟─e252d8ea-9763-4030-9440-59e3a63ec4c2
# ╠═590ea7ba-eeae-4f03-8189-a199714001fe
# ╟─5a2320c0-69c9-4ab2-b350-767b9db2d228
# ╠═d89b13de-0673-4cce-ab68-c3ae7bc89d5d
# ╟─b7a00261-686e-4c29-a11f-42764c60762d
# ╠═309eb417-5487-4362-990b-6221cdb8fc88
# ╠═97b63bdc-8655-4713-9908-d04ccc801932
# ╠═7299f93d-9318-4778-812f-1cb9246c2110
# ╠═64083d10-8846-4d18-8893-ff5d04a3941a
# ╟─0fcc8354-8139-4c72-acc1-c9378939e275
# ╠═24f2c0dc-d99b-4809-9080-103a48eb0083
# ╠═cf11045e-c8ef-48f0-8ab7-eb38b5119398
# ╟─b63c5e21-2218-42ba-80da-dcd74d0ddc38
# ╠═4e620b69-903b-4085-b45a-cbd8f12a5788
# ╠═59b01860-6c09-43d7-ae32-70561fec0e58
# ╠═4977c8c7-8a43-4ef0-ad18-2aaa9c147706
# ╠═998c5d5a-a28f-43ac-900f-508da56e9d3b
# ╠═09a701f0-ea44-419f-9de4-f5b4682e9d94
# ╠═8a653fd5-f5c5-4fc6-88d7-6eeb4da81216
# ╠═5a3f44fc-2016-49a9-9e32-30c9a2bff081
# ╠═67e55bab-f2cd-4181-b9a0-33b2181e207b
# ╠═d3256762-6488-4b11-8490-261a98fad0b2
# ╠═747e8dbb-d41a-442d-a483-bef308ef7f40
# ╠═e7dd85e6-e03e-48ad-96fd-6b9e7c4cb968
# ╠═6d64bbd1-dbe5-463d-9fd5-bf72fb92ad45
# ╠═235dfa07-9c56-47a2-a4a8-0b540625bc26
# ╠═941c7e7a-6b6d-48a8-8426-426b9ddcb9e1
# ╠═7ad2f7f3-e6b6-4707-ac81-22c51d4a8316
# ╠═6b705783-890d-49f2-b3a2-6228897e2fb4
# ╠═2a1daba5-9035-4a5c-b32d-c73c022135c6
# ╠═27b79f65-fd20-4970-927c-1fe5a453487c
# ╠═dbc9c3c2-e6c6-4d3c-81f0-3b8c3dcd4c02
# ╟─e6ed275d-93b4-4e9b-a05b-5d59915536b6
# ╠═97fa01d5-ecd4-453d-9963-9b7d50eb67cc
# ╠═38baae5d-3d15-4c79-840a-59635caf7636
# ╠═8915f5e6-e9cd-4ba4-a073-bdb1cb1f7890
# ╠═5bc565c6-25fd-44e1-8222-3741d0e13afb
# ╟─fb468dfe-4176-4b84-90a0-6346fe4c7475
# ╠═3b828c13-737a-4b5f-8f1f-02d4c7384b6a
# ╠═3b04bdf7-24d8-4269-8576-02ac2bad6c97
# ╠═88ef84f4-b1a6-4159-bbbf-f4f506859950
# ╠═b2b15ac7-1574-4999-bbc9-3536d3aa18b1
# ╠═09ca7aa2-70a7-4d1b-ae9f-bdd937bf5c94
# ╠═e0e8eb60-7ad7-4015-b226-56af6d405fad
# ╠═143d6069-6978-43ed-815f-5cf19de83fbf
# ╠═5433be8c-61fa-4262-839a-d9a137c05377
# ╠═c8d7966a-61c1-4ab2-856c-6ca5d10648cd
# ╠═358d297d-4def-44dc-b79d-6d23aa67d179
# ╠═c31ccfdc-3687-46cb-9eed-03791bde56d4
# ╠═75c4a277-9edb-4f47-bd5e-13dcb498bd96
# ╠═68b7380b-2514-4362-9504-552fb4e467f7
# ╠═f2a127b3-df43-486a-a34b-e289539a22bf
# ╠═f5806376-9ee8-4991-a84e-d650a149edd1
# ╠═7f0c719d-fbdf-455a-94f4-a1e899462c5d
# ╠═a8733344-6aba-47a5-bbc4-444fac026d40
