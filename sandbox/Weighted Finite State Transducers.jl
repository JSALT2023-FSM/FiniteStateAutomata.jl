### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ e2560be8-f3f9-11ed-0754-e3e1572c751d
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	
	# Register the FAST registry
	Pkg.Registry.add(
		RegistrySpec(url="https://gitlab.lisn.upsaclay.fr/fast/registry")
	)
	Pkg.add("Semirings")
	Pkg.add("PlutoUI")

	using LinearAlgebra
	using Semirings 
	using SparseArrays
	using PlutoUI
end

# ╔═╡ a8529ad0-6851-473f-bcc2-8bffd6ffc9a0
md"""
# Matrix-based representation of a WFST
*Lucas ONDEL YANG, LISN - CNRS, 2023*

Prototype of a linear algebra-based WFST library. 
"""

# ╔═╡ e7b33149-d195-4b9f-ae75-b520808e0928
TableOfContents()

# ╔═╡ 3e42b849-1a7e-41af-944d-6aed81965ba3
md"""
## Representing WFST 

Let ``\mathcal{A} = (Λ, Q, I, F, α, ω, E)`` be a WFST. Furthermore we assume that each arc ``e \in E`` is assigned to a unique id ``i \in \{1, \dots, |E|\}``. We encode the arcs with a (source) matrix ``\mathbf{S}`` and a (destination) matrix ``\mathbf{D}``. 

``\mathbf{S}`` is a ``|Q| \times |E|`` matrix where:
```math
	S_{ij} = \begin{cases}
		w[e] & \text{if there is an arc } e \text{ with id } j \text{ and leaving state } i \\
		\bar{0} & \text{otherwise}
	\end{cases}
```.

``\mathbf{D}`` is a ``|E| \times |Q|`` matrix where:
```math
	D_{ij} = \begin{cases}
		\bar{1} & \text{if there is an arc } e \text{ with id } i \text{ and going to state } j \\
		\bar{0} & \text{otherwise}
	\end{cases}
```.


The intial (resp. final) weights are encoded in a (sparse) vector ``\boldsymbol{\alpha}`` (resp. ``\boldsymbol{\omega}``).
""" 

# ╔═╡ 7d62a096-6f2a-40eb-914a-cb19fe998523
struct WFST{K,L<:Union{Int,Pair{Int,Int}}}
	α::AbstractVector{K}
	S::AbstractMatrix{K}
	D::AbstractMatrix{K}
	ω::AbstractVector{K}
	λ::Vector{L} # <- arc labels 
end

# ╔═╡ af8ccb1b-acdb-4b60-906e-d400e7f45a16
nstates(A::WFST) = size(A.S, 1)

# ╔═╡ eb0ac6e2-0ccf-4878-b037-a794ec430b82
narcs(A::WFST) = size(A.S, 2)

# ╔═╡ 946bac07-8939-45a4-9b13-94c49d3431e8
md"""
## Interfacing with OpenFST
"""

# ╔═╡ dbbd9409-c407-4403-9d1a-1a903623daa5
isyms = """
ϵ 0
a 1
b 2
c 3
"""

# ╔═╡ 34c98dad-3df6-444b-833e-7e7dd8e3e619
osyms = """
ϵ 0
x 1
y 2
z 3
"""

# ╔═╡ 377c3588-fa72-4815-a51c-50beba96b4e9
function symboltable(txt::AbstractString)
	symtable = Dict()
	for line in split(txt, "\n")
		tokens = split(line)
		length(tokens) == 0 && continue # skip emtpy lines
		label, id = tokens 
		symtable[parse(Int, id)] = label
	end
	symtable
end

# ╔═╡ 8131fd7c-ca50-490a-8525-fa6c863606dc
begin 
	function compile(wfst::AbstractString; 
					 semiring = LogSemiring{Float32},
					 acceptor = false)
		K = semiring
		arcslabel = acceptor ? Int[] : Pair{Int,Int}[]
		arcsweight = K[]
		arcssrc = Int[]
		arcsdest = Int[]
		fstates = Int[]
		fweights = K[]
		
		lines = split(wfst, "\n")
		for line in lines
			tokens = split(line)
			isempty(tokens) && continue 
			
			if 1 ≤ length(tokens) ≤ 2
				state = parse(Int, tokens[1]) + 1
				weight = length(tokens) == 1 ? one(K) : K(parse(Float64, tokens[2]))
				push!(fstates, state)
				push!(fweights, weight)
			else
				src = parse(Int, tokens[1]) + 1
				dest = parse(Int, tokens[2]) + 1
				push!(arcssrc, src)
				push!(arcsdest, dest)

				isym = parse(Int, tokens[3])
				if ! acceptor
					osym = parse(Int, tokens[4])
					weight = length(tokens) == 4 ? one(K) : K(parse(Float64, tokens[4]))
					push!(arcslabel, isym => osym)
				else
					weight = length(tokens) == 3 ? one(K) : K(parse(Float64, tokens[3]))
					push!(arcslabel, isym)
				end
				push!(arcsweight, weight)
			end
		end

		Qset = Set{Int}(arcssrc) ∪ Set{Int}(arcsdest)

		initstate = parse(Int, split(first(lines))[1]) + 1
		Q = length(Qset)
		A = length(arcslabel)
	
		WFST(
			sparsevec([initstate], one(K), Q),
			sparse(arcssrc, collect(1:A), arcsweight, Q, A),
			sparse(collect(1:A), arcsdest, one(K), A, Q),
			sparsevec(fstates, fweights, Q),
			arcslabel
		)
	end
end

# ╔═╡ 9b5251a7-c77e-4062-b274-b7c5ad966f97
compile(
	"""
	0 1 1 1 .5
	0 1 2 2 1.5
	1 2 3 3 2.5
	2 3.5
	"""
)

# ╔═╡ 53d3a8f7-d5c3-4363-816a-8846599e7101
compile(
	"""
	0 1 1 .5
	0 1 2 1.5
	1 2 3 2.5
	2 3.5
	""";
	acceptor = true
)

# ╔═╡ 20c14cc3-0264-42fe-be7b-399ed9debef1
compile(
	"""
	0 1 1 .5
	0 1 2 1.5
	1 2 3 2.5
	2 3.5
	1 12
	""";
	semiring = TropicalSemiring{Float32},
	acceptor = true
)

# ╔═╡ 3754c2ab-9134-4652-afb1-a024f37e4fa4
md"""
## Example 

### Creating a WFST
"""

# ╔═╡ 2b8f46d1-08a5-479a-8261-7dc6f662563a
K = ProbSemiring{Float32}

# ╔═╡ a659107e-df45-478b-962c-c43f1730b8d2
A = WFST(
	sparsevec([1], one(K), 3),
	sparse([1, 1, 2], [1, 2, 3], K[0.5, 1.5, 2.5], 3, 3),
	sparse([1, 2, 3], [2, 2, 3], one(K), 3, 3),
	sparsevec([3], K(3.5), 3),
	[1=>1, 2=>2, 3=>3]
)

# ╔═╡ e3848c7d-3c0d-4476-b948-17775377afdf
md"""
### Visualization
"""

# ╔═╡ 89eca3ce-e46a-4e57-a7ed-6a7ae849742c
sym1 = Dict(
	1 => "a",
	2 => "b",
	3 => "c"
)

# ╔═╡ d22ad598-d5ce-44fc-b9ee-f7299c7f1593
sym2 = Dict(
	1 => "x",
	2 => "y",
	3 => "z"
)

# ╔═╡ 93f41798-3fc0-4422-85c4-0c508c975849
(wfst=A, sym=sym1) # use "sym1" for both input and output labels

# ╔═╡ 05f2e1cd-3912-4b52-8004-5982e79d3ffe
(wfst=A, isym=sym1, osym=sym2)

# ╔═╡ 89afc39c-05fe-4452-a56c-d9317979702b
md"""
### Computing the WFST weight

Let be the matrix ``\mathbf{T} = \mathbf{S} \mathbf{D}``.

The weight of the WFST ``W(\mathcal{A})``, i.e. the ``\oplus``-sum of all its paths' weight, is given by
```math
W(\mathcal{A}) = \boldsymbol{\alpha}^\top (\mathbf{T}^0 + \mathbf{T}^1 + \mathbf{T}^2 + \dots) \boldsymbol{\omega}.
```

Which is trivially computed using dynamic programming:
```math
\begin{align}
	\mathbf{u}_0 &= \boldsymbol{\alpha} \\
	\mathbf{u}_n &= (\mathbf{u}_{n-1}^\top \mathbf{T})^\top \\
	W(\mathcal{A}) &= \mathbf{u}_0^\top \boldsymbol{\omega} \oplus \mathbf{u}_1^\top \boldsymbol{\omega} \oplus \mathbf{u}_2^\top \boldsymbol{\omega} \oplus \dots
\end{align}
```
"""

# ╔═╡ 425f9551-e38c-433e-a080-85038ff8d23c
function W(A::WFST)
	uₙ = A.α
	T = A.S * A.D
	
	retval = dot(uₙ, A.ω)
	n = 0
	while nnz(uₙ) > 0 
		n > nstates(A) && throw(ArgumentError("cyclic WFST"))
		n += 1

		# The next two lines of code are the bottleneck 
		# for MMI loss evaluation (and it's gradient).
		# Should be optimized at some point.
		uₙ = transpose(transpose(uₙ) * T)
		retval = retval ⊕ dot(uₙ, A.ω)
	end
	retval
end

# ╔═╡ 0125ad80-de6f-4d35-ab31-0a015ed894aa
W(A)

# ╔═╡ 7052e9db-6070-4f6e-8f98-24fc20bc175d
md"""
## Graphviz Backend for Visualization
"""

# ╔═╡ 0c84ee4e-1884-49ec-adc7-072647dd80d8
begin
	function _show_svg(io::IO, A::WFST; isym = Dict(), osym = Dict())
	    dotpath, dotfile = mktemp()
	    try
			fileio = IOContext(dotfile, :compact => true, :limit => true)
	        show(fileio, MIME("text/vnd.graphviz"), A; isym, osym)
	        close(dotfile)
	        svgpath, svgfile = mktemp()
	        try
	            run(`dot -Tsvg $(dotpath) -o $(svgpath)`)
	            xml = read(svgfile, String)
	            write(io, xml)
	        finally
	            close(svgfile)
	            rm(svgpath)
	        end
	    finally
	        rm(dotpath)
	    end
	end
	
	Base.show(io::IO, ::MIME"image/svg+xml", A::WFST) = _show_svg(io, A)
	
	Base.show(io::IO, ::MIME"image/svg+xml", tup::NamedTuple{(:wfst, :sym)}) = 
	_show_svg(io, tup.wfst; isym = tup.sym, osym = tup.sym)
	
	Base.show(io::IO, ::MIME"image/svg+xml", tup::NamedTuple{(:wfst, :isym, :osym)}) = 
	_show_svg(io, tup.wfst; isym = tup.isym, osym = tup.osym)
	
	function Base.show(io::IO, ::MIME"text/vnd.graphviz", A::WFST; isym = Dict(), osym = Dict())
		println(io, "Digraph {")
		println(io, "rankdir=LR;")
	
		Q = nstates(A)
		P = narcs(A)
		
		for q in 1:nstates(A)
			label = "$q"
			style = "solid"
			shape = "circle"
			
			if ! iszero(A.α[q])
				style = "bold"
				if ! isone(A.α[q])
					label = "$q/$(A.α[q])"
				end
			end

			if ! iszero(A.ω[q])
				shape = "doublecircle"
				if ! isone(A.ω[q])
					label = "$q/$(A.ω[q])"
				end
			end

			if ! iszero(A.α[q]) && ! isone(A.α[q]) && ! iszero(A.ω[q]) && ! isone(A.ω[q])
				println(io, q, " [label=\"$q/", A.α[q], "/", A.ω[q], "\", style=\"$style\", shape=\"$shape\"];")
			elseif ! iszero(A.α[q]) && ! isone(A.α[q])
				println(io, q, " [label=\"$q/", A.α[q], "\", style=\"$style\", shape=\"$shape\"];")
			elseif ! iszero(A.ω[q]) && ! isone(A.ω[q])
				println(io, q, " [label=\"$q/", A.ω[q], "\", style=\"$style\", shape=\"$shape\"];")
			else
				println(io, q, " [label=\"$q\", style=\"$style\", shape=\"$shape\"];")
			end
		end
	
		_get = (sym, l) -> get(sym, l, l)
		for (s, a, w, _, d) in zip(findnz(A.S)..., findnz(A.D)...)
			l = A.λ[a]
			label = l isa Pair ? join([_get(isym, first(l)), _get(osym, last(l))], ":") : _get(isym, first(l))
			println(io, s, "->", d, " [label=\"", label, "/", w, "\"];")
		end
	
		println(io, "}")
	end
end

# ╔═╡ 6e669195-4617-4344-874e-6a5aab490c93
show(stdout, MIME("text/vnd.graphviz"), A)

# ╔═╡ Cell order:
# ╟─a8529ad0-6851-473f-bcc2-8bffd6ffc9a0
# ╠═e2560be8-f3f9-11ed-0754-e3e1572c751d
# ╠═e7b33149-d195-4b9f-ae75-b520808e0928
# ╟─3e42b849-1a7e-41af-944d-6aed81965ba3
# ╠═7d62a096-6f2a-40eb-914a-cb19fe998523
# ╠═af8ccb1b-acdb-4b60-906e-d400e7f45a16
# ╠═eb0ac6e2-0ccf-4878-b037-a794ec430b82
# ╟─946bac07-8939-45a4-9b13-94c49d3431e8
# ╠═dbbd9409-c407-4403-9d1a-1a903623daa5
# ╠═34c98dad-3df6-444b-833e-7e7dd8e3e619
# ╠═377c3588-fa72-4815-a51c-50beba96b4e9
# ╠═8131fd7c-ca50-490a-8525-fa6c863606dc
# ╠═9b5251a7-c77e-4062-b274-b7c5ad966f97
# ╠═53d3a8f7-d5c3-4363-816a-8846599e7101
# ╠═20c14cc3-0264-42fe-be7b-399ed9debef1
# ╟─3754c2ab-9134-4652-afb1-a024f37e4fa4
# ╠═2b8f46d1-08a5-479a-8261-7dc6f662563a
# ╠═a659107e-df45-478b-962c-c43f1730b8d2
# ╟─e3848c7d-3c0d-4476-b948-17775377afdf
# ╠═6e669195-4617-4344-874e-6a5aab490c93
# ╠═89eca3ce-e46a-4e57-a7ed-6a7ae849742c
# ╠═d22ad598-d5ce-44fc-b9ee-f7299c7f1593
# ╠═93f41798-3fc0-4422-85c4-0c508c975849
# ╠═05f2e1cd-3912-4b52-8004-5982e79d3ffe
# ╟─89afc39c-05fe-4452-a56c-d9317979702b
# ╠═425f9551-e38c-433e-a080-85038ff8d23c
# ╠═0125ad80-de6f-4d35-ab31-0a015ed894aa
# ╟─7052e9db-6070-4f6e-8f98-24fc20bc175d
# ╠═0c84ee4e-1884-49ec-adc7-072647dd80d8
