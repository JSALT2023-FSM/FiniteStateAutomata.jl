### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 7911b2c4-183f-11ee-05e9-c9b311ee4780
begin
	using Pkg
	Pkg.develop(path="../../finitestateautomata.jl/")
	using Revise
	using SparseArrays
	using FiniteStateAutomata
    using PlutoUI, BenchmarkTools
	using ProfileCanvas
	using Profile
	using Random
end

# ╔═╡ 38eca186-37ed-4ae7-8341-3e28fee6ecc0
HTML("""
<!-- the wrapper span -->
<div>
	<button id="myrestart" href="#">Restart</button>
	
	<script>
		const div = currentScript.parentElement
		const button = div.querySelector("button#myrestart")
		const cell= div.closest('pluto-cell')
		console.log(button);
		button.onclick = function() { restart_nb() };
		function restart_nb() {
			console.log("Restarting Notebook");
		        cell._internal_pluto_actions.send(                    
		            "restart_process",
                            {},
                            {
                                notebook_id: editor_state.notebook.notebook_id,
                            }
                        )
		};
	</script>
</div>
""")

# ╔═╡ 2c3c6583-491c-4af9-a424-7d5f983ee982
symbols = Dict(1 => "a", 2 => "b", 3 => "c")

# ╔═╡ 438a3823-a65f-4ab3-9cf0-7ebf946c9226
S = TropicalSemiring{Float64}

# ╔═╡ 98d1b5bd-03e4-45be-bee3-1960f830253f
macro lencheck(l, vars...)
  exprs = Expr[]
  for var in vars
	varname = string(var)
	push!(exprs, :(
	  if length($(esc(var))) != $(esc(l))
		throw(DimensionError($varname, $(esc(l)), length($(esc(var)))))
  end
	))
  end
  Expr(:block, exprs...)
end

# ╔═╡ 13ddcb81-1c2a-4d2c-bc2e-0a7dcdb0d1fb
begin
	abstract type AbstractSparseMatrixCOO{Tv, Ti <: Integer} <: AbstractSparseMatrix{Tv, Ti} end
	
	mutable struct SparseMatrixCOO{Tv, Ti <: Integer} <: AbstractSparseMatrixCOO{Tv, Ti}
	  m::Int
	  n::Int
	  rows::Vector{Ti}
	  cols::Vector{Ti}
	  vals::Vector{Tv}
	
	  function SparseMatrixCOO{Tv, Ti}(
	    m::Integer,
	    n::Integer,
	    rows::Vector{Ti},
	    cols::Vector{Ti},
	    vals::Vector{Tv},
	  ) where {Tv, Ti <: Integer}
	    @noinline throwsz(str, lbl, k) =
	      throw(ArgumentError("number of $str ($lbl) must be ≥ 0, got $k"))
	    m < 0 && throwsz("rows", 'm', m)
	    n < 0 && throwsz("columns", 'n', n)
	    nnz = length(vals)
	    @lencheck nnz rows cols
	    new(Int(m), Int(n), rows, cols, vals)
	  end
	end

	Base.size(A::SparseMatrixCOO) = (getfield(A, :m), getfield(A, :n))

	# function Base.show(io::IO, ::MIME"text/plain", S::AbstractSparseMatrixCOO)
	#   	xnnz = nnz(S)
	#   m, n = size(S)
	#   print(io, m, "×", n, " ", typeof(S), " with ", xnnz, " stored ", xnnz == 1 ? "entry" : "entries")
	#   if xnnz != 0
	#     print(io, ":\n")
	#     show(IOContext(io, :typeinfo => eltype(S)), S)
	#   end
	# end

	function Base.print_array(io::IO, S::AbstractSparseMatrixCOO)
	    if max(size(S)...) < 100
			D = todense(S)
	        Base.print_matrix(io, D)
	    else
	        _show_with_braille_patterns(io, S)
	    end
	end

	function todense(A::SparseMatrixCOO) 
		D = zeros(eltype(A),(A.m,A.n))
	  	for k = 1:nnz(A)
	    	i, j, v = A.rows[k], A.cols[k], A.vals[k]
	    	D[i,j] = v
	    end
		D
	end

	nnz(A::SparseMatrixCOO) = length(A.vals)
	nnz(A::Nothing) = 0

	function Base.getindex(A::AbstractSparseMatrixCOO{Tv, Ti}, i0::Integer, i1::Integer) where {Tv, Ti}
	  m, n = size(A)
	  (1 ≤ i0 ≤ m && 1 ≤ i1 ≤ n) || throw(BoundsError())
	  for k = 1:nnz(A)
	    i, j, v = A.rows[k], A.cols[k], A.vals[k]
	    if i == i0 && j == i1
	      return v
	    end
	  end
	end

	function csc2coo(A)
		colptr = A.colptr
		rowval = A.rowval
		vals = A.nzval
		rows = Vector{Int}()
		cols = Vector{Int}()
		for i in 1:length(colptr)-1
			for j in colptr[i]:colptr[i+1]-1
				push!(rows,rowval[j])
				push!(cols,i)
			end
		end
		SparseMatrixCOO{eltype(vals),Int}(A.m,A.n,rows, cols, vals)
	end
end

# ╔═╡ 26b4c706-09ec-463e-b9b5-37e473b17461
begin
	function KronCOO(A,B)
		rows = Vector{Int}()
		cols = Vector{Int}()
		S = eltype(A)
		vals = Vector{S}()
		ma, mb, na, nb = A.m, B.m, A.n, B.n 
		for (i,j,a) in zip(A.rows, A.cols, A.vals)
			for (k,l,b) in zip(B.rows, B.cols, B.vals)
				c = a*b
				push!(rows,(i-1)*mb+k)
				push!(cols,(j-1)*nb+l)
				push!(vals, c)
			end
		end
		SparseMatrixCOO{S, Int}(ma*mb, na*nb, rows, cols, vals)
	end
	
	Base.kron(A::AbstractSparseMatrixCOO{S}, B::AbstractSparseMatrixCOO{S}) where S =
	    KronCOO(A, B)

	# using matrix multiplication with csc, returns coo
	# missing matmul in csr from semimodules 
	function Base.:*(A::AbstractSparseMatrixCOO, B::AbstractSparseMatrixCOO)
		csc2coo(sparse(A.rows,A.cols,A.vals, A.m, A.n)*sparse(B.rows,B.cols,B.vals, B.m, B.n))
	end
end

# ╔═╡ 2c3ae5ef-1e9a-420e-973f-fa76b4ba2734
md"Testint COO"

# ╔═╡ 02954571-531a-455a-ace9-050cd42c7d89
SparseMatrixCOO{Float64, Int64}(8,8,[1, 1, 2, 3],[1, 7, 2, 4],[1.0, 2.0, -5.0, 3.0])

# ╔═╡ 6275c609-0cc0-4c4e-8fba-59442299fcab
md"Two level 2D sparse COO, kron on the outer level, mat mul in the inner level"

# ╔═╡ 6bef7c64-2c87-4bdf-80b1-148ac46c479c
begin
	Random.seed!(1234)
	nv = 4
	ns = 3
	ms = 5
	mval = 8
	V = Vector{SparseMatrixCOO{S, Int}}()
	for i in 1:mval
		vals = [S(v) for v in rand(nv)]
		Stmp = SparseMatrixCOO{S, Int}(nv,nv,rand(1:ns,nv),rand(1:ns,nv),vals)
		push!(V,Stmp)
	end
	M = SparseMatrixCOO{SparseMatrixCOO{S, Int}, Int}(ms,ms,rand(1:ms,mval),rand(1:ms,mval),V);
	k = kron(M, M);
	"Missing show function for 2D coo"
end

# ╔═╡ b64e8b2b-9c69-4cce-9595-a04265c9a299
md"
## FST Composition Example
"

# ╔═╡ fd4f5379-7d6b-4962-8f07-ddedacf565a6
begin
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
	B = VectorFST(
		[
			Arc{S}[(2, 2, 3, S(0.3))],
			Arc{S}[(3, 1, 2, S(0.4))],
			Arc{S}[(3, 1, 2, S(0.6))],
		],
		1,
		S[zero(S),zero(S),S(0.7)]
	)
end

# ╔═╡ ddff62e9-b7ed-43b3-beaf-c46c8ebe61ec
draw(A; symbols)

# ╔═╡ 2d6277aa-d73f-43be-9a91-2e82218a5587
draw(B; symbols)

# ╔═╡ bfb7712b-59cb-4666-bbbd-406c7ee8a5e9
# helper functions TODO: use less
begin
	function coo2arcs(A)
		t = Vector{Tuple{Int,Int,Int,S}}
		nstates = A.m
		states = Vector{t}()
		for i in 1:nstates
			push!(states,Vector{t}())
		end
	
		for (s,d,m) in zip(A.rows, A.cols, A.vals)
			for (i,j,v) in zip(m.rows, m.cols, m.vals)
				push!(states[s],(d,i,j,v))
			end
		end
		states
	end
	function vector2dict(A)
		# VectorFST to DictFST
		Q = numstates(A)
		d = Dict{Int,Dict}()
		for (s,arcs) in enumerate(A.arcs)
			for arc in arcs
				if !haskey(d,s)
					d[s] = Dict{Int,Vector{Tuple{Int,Int,S}}}()
				end
				if !haskey(d[s], arc[1])
					d[s][arc[1]] = Vector{Tuple{Int,Int,S}}()
				end
				push!(d[s][arc[1]], arc[2:4])
			end
		end
		d
	end
	function dict2coo(dict_fst, nstates, nsyms)
		state_rows = Vector{Int}()
		state_cols = Vector{Int}()
		state_vals = Vector{SparseMatrixCOO{S, Int}}()
		for (s, sdic) in dict_fst
		    for (d, ddic) in sdic
				push!(state_rows,s)
				push!(state_cols,d)
				label_rows = Vector{Int}()
				label_cols = Vector{Int}()
				label_vals = Vector{S}()
				for (i,j,v) in ddic
					push!(label_rows, i)
					push!(label_cols, j)
					push!(label_vals, v)
				end
				push!(state_vals, SparseMatrixCOO{S, Int}(nsyms, nsyms, label_rows, label_cols, label_vals))
			end
		end
		SparseMatrixCOO{SparseMatrixCOO{S, Int}, Int}(nstates, nstates, state_rows, state_cols, state_vals)
	end;
end

# ╔═╡ 9ec6ce15-33c2-4155-9572-8097d4caac58
function sparse_composition_matmutl(A,B, nsymbols)
	cooA = dict2coo(vector2dict(A), numstates(A),  nsymbols)
	cooB = dict2coo(vector2dict(B), numstates(B),  nsymbols)
	cooC = kron(cooA, cooB)
	arcsC = coo2arcs(cooC)
	# TODO initial state
	VectorFST(arcsC, 1 , kron(A.finalweights, B.finalweights) )
end

# ╔═╡ 1f665bf9-99a2-4244-a0a8-87d5c8a1fcb0
begin
	C = sparse_composition_matmutl(A,B, length(symbols))
	draw(C; symbols)
end

# ╔═╡ Cell order:
# ╟─38eca186-37ed-4ae7-8341-3e28fee6ecc0
# ╠═7911b2c4-183f-11ee-05e9-c9b311ee4780
# ╠═2c3c6583-491c-4af9-a424-7d5f983ee982
# ╠═438a3823-a65f-4ab3-9cf0-7ebf946c9226
# ╠═98d1b5bd-03e4-45be-bee3-1960f830253f
# ╠═13ddcb81-1c2a-4d2c-bc2e-0a7dcdb0d1fb
# ╠═26b4c706-09ec-463e-b9b5-37e473b17461
# ╟─2c3ae5ef-1e9a-420e-973f-fa76b4ba2734
# ╠═02954571-531a-455a-ace9-050cd42c7d89
# ╟─6275c609-0cc0-4c4e-8fba-59442299fcab
# ╠═6bef7c64-2c87-4bdf-80b1-148ac46c479c
# ╟─b64e8b2b-9c69-4cce-9595-a04265c9a299
# ╠═fd4f5379-7d6b-4962-8f07-ddedacf565a6
# ╠═ddff62e9-b7ed-43b3-beaf-c46c8ebe61ec
# ╠═2d6277aa-d73f-43be-9a91-2e82218a5587
# ╠═bfb7712b-59cb-4666-bbbd-406c7ee8a5e9
# ╠═9ec6ce15-33c2-4155-9572-8097d4caac58
# ╠═1f665bf9-99a2-4244-a0a8-87d5c8a1fcb0
