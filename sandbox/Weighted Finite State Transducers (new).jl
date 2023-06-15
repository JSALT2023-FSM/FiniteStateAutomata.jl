### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ a2acb152-0475-11ee-11af-8b864963e161
begin
	using Revise
	using Pkg
	Pkg.develop(path="..")
	Pkg.add("PlutoUI")
	using FiniteStateAutomata
	using Semirings
	using SparseArrays
	using PlutoUI
end

# ╔═╡ 69efb345-afe4-435d-afd4-30c94109a7aa
TableOfContents()

# ╔═╡ 7d604912-13ca-4882-8420-6b7b681a95e7
TableOfContents(title="📚 Table of Contents", indent=true, depth=4)

# ╔═╡ 850192a5-9578-4e2c-abb3-b8fa46cfda72
md"""
## Creating WFST
"""

# ╔═╡ bc5265f3-31d1-4387-8147-9281a3cd1da9
syms1 = symboltable(
	"""
	ϵ 0
	a 1
	b 2 
	c 3
	d 4
	e 5
	"""
)

# ╔═╡ ce17161c-3315-4be8-aaca-af7b36b18175
syms2 = symboltable(
	"""
	ϵ 0
	q 1
	r 2
	s 3
	"""
)

# ╔═╡ ff026ed9-9062-420a-9b31-b026031865b1
syms3 = symboltable(
	"""
	ϵ 0
	f 1
	g 2
	h 3
	j 4
	"""
)

# ╔═╡ 16f91e50-6bcd-4ce3-9d43-a91d2e6396ee
fst = compile(
	"""
	0 1 1 1 .5
	0 1 2 2 1.5
	1 2 3 3 2.5
	2 3.5
	""";
	semiring = TropicalSemiring{Float32}, # Default LogSemiring{Float32}
	openfst_compat = true
)

# ╔═╡ ce84a5d1-4ffe-4e9a-993d-e6c97f5d8b59
compile(
	"""
	1 2 1 1 .5
	1 2 2 2 1.5
	2 3 3 3 2.5
	3 3.5
	""";
	semiring = TropicalSemiring{Float32}, 
)

# ╔═╡ 8e1cca47-0677-4d13-a080-08c64a09d14e
fsa = compile(
	"""
	1 2 1 .5
	1 2 2 1.5
	2 3 3 2.5
	3 3.5
	""";
	semiring = TropicalSemiring{Float32}, # Default LogSemiring{Float32}
	acceptor = true
)

# ╔═╡ 9bc55c91-fe21-41ef-bee7-6f67621fe07a
md"""
## Exporting FST
"""

# ╔═╡ fc8f3443-f898-4c0f-aea3-f0cfda5af600
print(fst)

# ╔═╡ 54b8665e-ff95-4fce-89b9-106296d063f8
print(fst; openfst_compat = true)

# ╔═╡ d80caeac-40c6-4e3b-b0c8-bf3c63f60cf5
print(fsa)

# ╔═╡ ac007b42-991b-42c2-8395-7cba132c8abc
md"""
### Visualization
"""

# ╔═╡ 61cabe58-2596-4c68-977b-e1a94528a44b
draw(fst)

# ╔═╡ 9e9551e6-f0ef-4bdc-a42f-c3570f9e92ea
draw(fst; isymbols = syms1, osymbols = syms2)

# ╔═╡ 26761703-7a3e-460c-a359-230a16b33906
fst |> summary

# ╔═╡ ddf8b433-10ee-4783-be3f-ad716b252c22
fsa |> summary

# ╔═╡ 583d264d-e498-4a19-952b-fb50bebbf516
md"""
## Basic operations

### Arc Filtering
"""

# ╔═╡ 28187f82-f099-4b7a-b0f4-71a6f3fc7daa
md"""
### State Filtering 
"""

# ╔═╡ 94cf6cf1-8c1b-4bd4-83c6-a4a1fce12854
md"""
### Project 
"""

# ╔═╡ 7d1c8038-ea43-4343-a992-b72cde854251
draw(fst; isymbols = syms1, osymbols = syms2)

# ╔═╡ a94fc98e-e487-412e-b45d-086eaf863004
draw(project(fst); isymbols = syms1, osymbols = syms1)

# ╔═╡ 6a5d2bf7-033d-4790-8938-7e8829ab1c18
draw(project(fst; project_output = true); isymbols = syms2, osymbols = syms2)

# ╔═╡ f7401134-9a29-4d1f-ba82-53773d0b5986
md"""
### Reverse
"""

# ╔═╡ f7d117ff-1213-4656-b5b5-2f603c89bf11
draw(fst; isymbols = syms1, osymbols = syms2)

# ╔═╡ 5e3faa50-ed40-49d4-bb26-15c1287ffe90
draw(reverse(fst); isymbols = syms1, osymbols = syms2)

# ╔═╡ de319274-ebbf-4b2b-9895-f1c901d0a21b
md"""
### Relabeling
"""

# ╔═╡ 0996e65c-8e21-4fd3-b961-2233fe3791d0
draw(fst)

# ╔═╡ 9bad9673-af64-4e44-b500-5a6b1fa9a1ad
rfst = relabel(fst) do l
	first(l) => last(l) + 1
end 

# ╔═╡ 186fb667-1ea6-499d-94b5-b9a4034d2dea
draw(rfst)

# ╔═╡ 15884f32-ac34-41fe-8335-a54a708f6216
md"""
### Inversion
"""

# ╔═╡ 45f4f380-b60e-428f-ab1c-7646b9960d10
draw(rfst)

# ╔═╡ b22fbd98-3b80-4dc7-8fdd-a4c7fdeb6a9e
draw(inv(rfst))

# ╔═╡ e622fe60-32d0-4433-ad5d-e5c5edabda96
md"""
### Tensor Product of FST
"""

# ╔═╡ 7ed60842-5f64-4415-a36c-33739203c2d4
md"""
### Composition 

#### Example 1: ϵ-free
"""

# ╔═╡ 409e679e-33e6-4b2c-9866-97dbb435d097
begin
	A = compile(
		"""
		1 2 1 1 1
	 	2 2 3 3 1
		1 3 1 2 2.5
		2
		3 2.5
		""";
		semiring = TropicalSemiring{Float32}
	)
	
	draw(A; isymbols = syms1, osymbols = syms2, openfst_compat = true)
end

# ╔═╡ 29d90b34-75fe-4138-8ea9-80f2dcde6ae1
fA1 = filterarcs(A) do (s, d, l, w)
	first(l) < 3
end

# ╔═╡ 152b52b3-11bd-436d-9d3d-4eecb6bc7072
draw(fA1; isymbols = syms1, osymbols = syms2)

# ╔═╡ f0a76cfe-e4eb-4e98-8b35-e2f257e68dea
fA2 = filterstates(A) do (q, iw, fw)
	q < 3
end 

# ╔═╡ 4fc4a673-4f93-47ed-95d2-37f316b1898e
draw(fA2; isymbols = syms1, osymbols = syms2)

# ╔═╡ 74862e5b-f860-4458-8e69-d0917e31cfe8
begin
	B = compile(
		"""
		1 2 1 1 1
		1 3 2 3 3
		2 3 3 2 2.5
		3 3 3 4 1.5
		3 2
		""";
		semiring = TropicalSemiring{Float32}
	)
	
	draw(B; isymbols = syms2, osymbols = syms3, openfst_compat = true)
end

# ╔═╡ f20c8851-7384-4a7c-a570-b8796e23f07c
draw(kron(A, B))

# ╔═╡ af3c6894-4a3c-4c97-8850-125d8c98e998
draw(A ∘ B; isymbols = syms1, osymbols = syms3)

# ╔═╡ 7d75c798-75d8-41e3-b055-b2e55abca9cc
md"""
#### Example 2: input and output ϵ labels
"""

# ╔═╡ 5943412c-b4aa-4b3c-b664-04e218927ae0
begin
	X = compile(
		"""
		1 2 1 1
		2 3 2 0
		3 4 3 0
		4 5 4 4 
		5
		"""
	)
	draw(X; isymbols = syms1, osymbols = syms1)
end

# ╔═╡ d9232933-311e-4607-98ca-683d5c145da0
begin
	Y = compile(
		"""
		1 2 1 4
		2 3 0 5
		3 4 4 1
		4 
		"""
	)
	draw(Y; isymbols = syms1, osymbols = syms1)
end	

# ╔═╡ 4356e0c8-3838-40c5-add9-7095e379718f
draw(X ∘ Y; isymbols = syms1, osymbols = syms1)

# ╔═╡ 41866c51-e609-4a8e-aaf5-86058e724e1c
C = kron(X, Y)

# ╔═╡ 5ccf3a8d-5cf6-4889-8f9c-2731b3dbd546
filterarcs(C) do (s, d, l, w)
        lA, lB = l
        last(lA) == first(lB)
end |> connect |> draw

# ╔═╡ 11bdc86e-d644-4995-88b9-cb08bbbf4688
party_fst = compile(
	"""
	1 2 1 1 
	1 2 0 1
	2 3 3 2 
	2 4 2 2 
	3 5 3 3 
	3 4 3 0 
	4 4 2 0
	4 4 3 0
	4 6 0 4
	5 
	6
	""";
) 

# ╔═╡ c8298bad-83ca-46a1-af87-ae76cf424b08
draw(
	party_fst; 
	isymbols = Dict(0 => "ϵ", 1 => "🍟",  2 => "🍺", 3 => "🍷"),
	osymbols = Dict(0 => "ϵ", 1 => "😐",  2 => "🙂", 3 => "😴", 4 => "🤮"),
)

# ╔═╡ 7b9e8e66-5191-426c-af7e-6cfd8fc84266
beer_fst = compile(
	"""
	1 2 1 1 0.5
	1 2 0 0 0.5
	2 3 2 2 0.3
	3 3 2 2 0.8
	2 4 3 3 0.7
	4 4 3 3 0.6
	3 0.2 
	4 0.4
	""";
) 

# ╔═╡ d6863b2e-647f-4345-b4e9-1dcd6f85290d
draw(
	beer_fst; 
	isymbols = Dict(0 => "ϵ", 1 => "🍟",  2 => "🍺", 3 => "🍷"),
	osymbols = Dict(0 => "ϵ", 1 => "🍟",  2 => "🍺", 3 => "🍷"),
)

# ╔═╡ 39167e85-fe2c-40e5-b86a-d3e22278bfc8
draw(
	beer_fst ∘ party_fst; 
	isymbols = Dict(0 => "ϵ", 1 => "🍟", 2 => "🍺", 3 => "🍷"),
	osymbols = Dict(0 => "ϵ", 1 => "😐",  2 => "🙂", 3 => "😴", 4 => "🤮"),
)

# ╔═╡ Cell order:
# ╠═a2acb152-0475-11ee-11af-8b864963e161
# ╠═69efb345-afe4-435d-afd4-30c94109a7aa
# ╠═7d604912-13ca-4882-8420-6b7b681a95e7
# ╟─850192a5-9578-4e2c-abb3-b8fa46cfda72
# ╠═bc5265f3-31d1-4387-8147-9281a3cd1da9
# ╠═ce17161c-3315-4be8-aaca-af7b36b18175
# ╠═ff026ed9-9062-420a-9b31-b026031865b1
# ╠═16f91e50-6bcd-4ce3-9d43-a91d2e6396ee
# ╠═ce84a5d1-4ffe-4e9a-993d-e6c97f5d8b59
# ╠═8e1cca47-0677-4d13-a080-08c64a09d14e
# ╟─9bc55c91-fe21-41ef-bee7-6f67621fe07a
# ╠═fc8f3443-f898-4c0f-aea3-f0cfda5af600
# ╠═54b8665e-ff95-4fce-89b9-106296d063f8
# ╠═d80caeac-40c6-4e3b-b0c8-bf3c63f60cf5
# ╟─ac007b42-991b-42c2-8395-7cba132c8abc
# ╠═61cabe58-2596-4c68-977b-e1a94528a44b
# ╠═9e9551e6-f0ef-4bdc-a42f-c3570f9e92ea
# ╠═26761703-7a3e-460c-a359-230a16b33906
# ╠═ddf8b433-10ee-4783-be3f-ad716b252c22
# ╟─583d264d-e498-4a19-952b-fb50bebbf516
# ╠═29d90b34-75fe-4138-8ea9-80f2dcde6ae1
# ╠═152b52b3-11bd-436d-9d3d-4eecb6bc7072
# ╟─28187f82-f099-4b7a-b0f4-71a6f3fc7daa
# ╠═f0a76cfe-e4eb-4e98-8b35-e2f257e68dea
# ╠═4fc4a673-4f93-47ed-95d2-37f316b1898e
# ╟─94cf6cf1-8c1b-4bd4-83c6-a4a1fce12854
# ╠═7d1c8038-ea43-4343-a992-b72cde854251
# ╠═a94fc98e-e487-412e-b45d-086eaf863004
# ╠═6a5d2bf7-033d-4790-8938-7e8829ab1c18
# ╟─f7401134-9a29-4d1f-ba82-53773d0b5986
# ╠═f7d117ff-1213-4656-b5b5-2f603c89bf11
# ╠═5e3faa50-ed40-49d4-bb26-15c1287ffe90
# ╟─de319274-ebbf-4b2b-9895-f1c901d0a21b
# ╠═0996e65c-8e21-4fd3-b961-2233fe3791d0
# ╠═9bad9673-af64-4e44-b500-5a6b1fa9a1ad
# ╠═186fb667-1ea6-499d-94b5-b9a4034d2dea
# ╟─15884f32-ac34-41fe-8335-a54a708f6216
# ╠═45f4f380-b60e-428f-ab1c-7646b9960d10
# ╠═b22fbd98-3b80-4dc7-8fdd-a4c7fdeb6a9e
# ╟─e622fe60-32d0-4433-ad5d-e5c5edabda96
# ╠═f20c8851-7384-4a7c-a570-b8796e23f07c
# ╟─7ed60842-5f64-4415-a36c-33739203c2d4
# ╟─409e679e-33e6-4b2c-9866-97dbb435d097
# ╟─74862e5b-f860-4458-8e69-d0917e31cfe8
# ╠═af3c6894-4a3c-4c97-8850-125d8c98e998
# ╟─7d75c798-75d8-41e3-b055-b2e55abca9cc
# ╟─5943412c-b4aa-4b3c-b664-04e218927ae0
# ╟─d9232933-311e-4607-98ca-683d5c145da0
# ╠═4356e0c8-3838-40c5-add9-7095e379718f
# ╠═41866c51-e609-4a8e-aaf5-86058e724e1c
# ╠═5ccf3a8d-5cf6-4889-8f9c-2731b3dbd546
# ╠═11bdc86e-d644-4995-88b9-cb08bbbf4688
# ╠═c8298bad-83ca-46a1-af87-ae76cf424b08
# ╟─7b9e8e66-5191-426c-af7e-6cfd8fc84266
# ╠═d6863b2e-647f-4345-b4e9-1dcd6f85290d
# ╠═39167e85-fe2c-40e5-b86a-d3e22278bfc8
