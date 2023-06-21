### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ ebb8cdee-0d15-11ee-0ff7-c9298a64b6ca
begin
	using Pkg
	Pkg.develop(path="../../finitestateautomata.jl")
	Pkg.develop(path="../../sparsesemimodules.jl")

	Pkg.add([
		PackageSpec(name="LogExpFunctions"),
		PackageSpec(url="https://gitlab.lisn.upsaclay.fr/fast/semirings.jl" , rev="jsalt2023"),
		PackageSpec(name="PlutoUI"),
		PackageSpec("CUDA")
	])

	using Revise
	using LinearAlgebra
	using FiniteStateAutomata
	using Semirings
	using PlutoUI
	using SparseSemimodules

	TableOfContents()
end

# ╔═╡ 9d6e8b9c-1e12-48ab-9043-86d24779ed8f
S = LogSemiring{Float32,1}

# ╔═╡ 0c7aef90-6533-4838-9373-ab07cc4ebfb0
A = SparseFST(
	SparseMatrices(
		sparse([1], [2], S(0.5), 3, 3),
		sparse([1], [2], S(1.5), 3, 3),
		sparse([2], [3], S(2.5), 3, 3),
	),
	sparsevec([1], one(S), 3),
	sparsevec([3], S(3.5), 3),
	[1 => 1, 2 => 2, 3 => 3]
)

# ╔═╡ 61c3400f-9a05-4462-8974-c4585817dea5
semiring(A)

# ╔═╡ a1b3b559-0862-4e18-80f5-c3323cbc4490
nstates(A)

# ╔═╡ fc7b5088-c0a2-441f-8fc0-dd8783e68004
arcs(A)

# ╔═╡ 644a1edd-9c75-4673-952a-50a6cfadc87a
states(A)

# ╔═╡ 211b42e3-46b0-49e4-a5ed-eb242b9b9709
draw(A; symbols = Dict(1 => "a", 2 => "b", 3 => "c"))

# ╔═╡ 7008d26d-d8cf-4366-83d9-cdbdc21b3a1d
B = compile(
	"""
	0 1 1 1 .5
	0 1 2 2 1.5
	1 2 3 3 2.5
	2 3.5
	""";
	openfst_compat = true
)

# ╔═╡ 2fd2500d-5c3c-4c3d-8cd4-b37f51bfd248
draw(B; symbols = Dict(1 => "a", 2 => "b", 3 => "c"))

# ╔═╡ 02113196-dd49-4580-b2ed-108c634cf980
md"""
## Direct sum

Let ``\mathbf{A}, \mathbf{B} \in S^{Q \times Q \times P}``. The direct sum ``\bar{\oplus}`` is given by:
```math
\begin{align}
	\mathbf{C} &= \mathbf{A} \bar{\oplus} \mathbf{B} &\implies
	\mathbf{C}_i &= \begin{bmatrix}
		\mathbf{A}_i & \mathbf{0} \\
		\mathbf{0} & \mathbf{B}_i
	\end{bmatrix} & \forall i \in \{1, \dots, P \}
\end{align}
```
"""

# ╔═╡ 3e9e7759-0747-4e13-bc14-b03c6ea94723
md"""
## Tensor product 

Let ``\mathbf{A}, \mathbf{B} \in S^{Q \times Q \times P}``. The tensor product is defined  ``\bar{\otimes}`` is given by:
```math
\begin{align}
	\mathbf{C} &= \mathbf{A} \bar{\otimes} \mathbf{B} &\implies
	\mathbf{C}_i &= \mathbf{A}_i \otimes \mathbf{B}_i = \begin{bmatrix}
		a_{i,11} \mathbf{B}_i & a_{i,12} \mathbf{B}_i & \dots \\
		a_{i,21} \mathbf{B}_i & \ddots \\
		\vdots
	\end{bmatrix} & \forall i \in \{1, \dots, P \}
\end{align}
```
"""

# ╔═╡ Cell order:
# ╠═ebb8cdee-0d15-11ee-0ff7-c9298a64b6ca
# ╠═9d6e8b9c-1e12-48ab-9043-86d24779ed8f
# ╠═0c7aef90-6533-4838-9373-ab07cc4ebfb0
# ╠═61c3400f-9a05-4462-8974-c4585817dea5
# ╠═a1b3b559-0862-4e18-80f5-c3323cbc4490
# ╠═fc7b5088-c0a2-441f-8fc0-dd8783e68004
# ╠═644a1edd-9c75-4673-952a-50a6cfadc87a
# ╠═211b42e3-46b0-49e4-a5ed-eb242b9b9709
# ╠═7008d26d-d8cf-4366-83d9-cdbdc21b3a1d
# ╠═2fd2500d-5c3c-4c3d-8cd4-b37f51bfd248
# ╟─02113196-dd49-4580-b2ed-108c634cf980
# ╟─3e9e7759-0747-4e13-bc14-b03c6ea94723
