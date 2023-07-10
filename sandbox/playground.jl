### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ fdf50f2e-1ef7-11ee-0d03-db187df6536b
begin
	using Pkg
	Pkg.develop(path="../../finitestateautomata.jl/")	
	using Revise
	using FiniteStateAutomata
    using PlutoUI, BenchmarkTools
	using Profile, ProfileCanvas
	using NPZ
end

# ╔═╡ dca59cae-aa9e-4dce-92b5-c4a2cca199ac
X = zeros(10,10,10)

# ╔═╡ f464e7ab-a831-473f-ad07-7f2de644956c
CartesianIndices((10,10))

# ╔═╡ 8fadec26-cb9e-4533-ad8f-fffb376fd170
2 == 2 == 2

# ╔═╡ 46a6da30-8b40-4b94-a7a1-e01a428813cf
for c in CartesianIndices((10,10))
	i, j = c[1],c[2]
	print(i,j)
end

# ╔═╡ b5f54187-2e6e-49f9-b327-6bfe501708c7
S[0.1,0.2,0.3]*S[0.1,0.2,0.3]

# ╔═╡ a239c43c-9c10-406f-8f7b-c1b7d383dd3a
symbols = Dict(1 => "a", 2 => "b", 3 => "c")

# ╔═╡ 145edb8f-9ffa-4e95-ba48-6c7dcdfeecd1
nsymbols = length(symbols)

# ╔═╡ efc53914-1c4e-4195-ac63-553f020c5f9b
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

# ╔═╡ db03e5e8-5e31-4f7d-88f6-e3096f56b5c4
function kron_coo2dict(A,B)  
    S = eltype(A)
    D = Dict()
    ma, mb, na, nb = A.m, B.m, A.n, B.n 
    for (i,j,a) in zip(A.rows, A.cols, A.vals)
        for (k,l,b) in zip(B.rows, B.cols, B.vals)
            c = a*b
			D[((i-1)*mb+k,(j-1)*nb+l )] = c
        end
    end
	D
end

# ╔═╡ 2513f1a1-92cd-44cb-afe9-a6ef5dddffb1
function sum_dicts(A,B)
	D = Dict()
	for a in keys(A)
		if !haskey(D,a)
			D[a]=0
		end
		D[a]+=A[a]
	end
	for b in keys(B)
		if !haskey(D,a)
			D[a]=0
		end
		D[a]+=A[a]
	end
end

# ╔═╡ 5c2e4fda-4304-4e6d-b98a-b51053f4b8ef
function my_sparse_composition_lod(A, B, nsymbols)	
	S = semiring(A)
    M = SparseMatrixCOO{S,Int}	
	
	cooA = dict2coo(vector2dict_lod(A), nsymbols, numstates(A), S)
	cooB = dict2coo(vector2dict_lod(B), nsymbols, numstates(B), S)
	
	Q = numstates(A) * numstates(B)
	
    label_rows = Vector{Int}()
    label_cols = Vector{Int}()
    label_vals = Vector{M}()

	kresults_dict = Dict()
	Threads.@threads for c in CartesianIndices((nsymbols,nsymbols,nsymbols))
		x,y,z = c[1],c[2],c[3]
		if hasitem(cooA,x,y) && hasitem(cooB,y,z)					
			kronresult = kron_coo2dict(cooA[x,y], cooB[y,z])
			if !haskey(kresults_dict, (x,z))
				kresults_dict[(x,z)] = []
			end
			push!(kresults_dict[(x,z)], kronresult)				
		end
	end

	result = Dict()
	for (label_pair, list_of_results) in kresults_dict
		
		buffer = Dict()		
		for kronresult in list_of_results
			for (newstates, value) in kronresult
				if !haskey(buffer, newstates)
					buffer[newstates]=S(0)
				end
				buffer[newstates]+=value
			end
		end
		result[label_pair] = buffer
	end
	result
	# cooC = SparseMatrixCOO{M, Int}(nsymbols, nsymbols, label_rows, label_cols, label_vals)
	# arcsC = coo_lod2arcs(cooC, Q, S)

	# initialA = zeros(S,numstates(A))
	# initialB = zeros(S,numstates(B))
	# initialA[A.initstate] = S(1)
	# initialB[B.initstate] = S(1)
	# initial = findfirst(x->x!=zero(S),kron(initialA, initialB))
	# final = kron(A.finalweights, B.finalweights) 

	# VectorFST(arcsC, initial, final )
	
end

# ╔═╡ 4b777a57-8437-4486-b6e4-184c6801f1ae
my_sparse_composition_lod(A,B, nsymbols)	


# ╔═╡ 3e4c65ed-8ab8-4d88-8fef-24270f9946a7
function my2_sparse_composition_lod(cooA, cooB, S, Q, nsymbols)	
    M = SparseMatrixCOO{S,Int}	
	
    label_rows = Vector{Int}()
    label_cols = Vector{Int}()
    label_vals = Vector{M}()

	kresults_dict = Dict()
	Threads.@threads for c in CartesianIndices((nsymbols,nsymbols,nsymbols))
		x,y,z = c[1],c[2],c[3]
		if hasitem(cooA,x,y) && hasitem(cooB,y,z)					
			kronresult = kron_coo2dict(cooA[x,y], cooB[y,z])
			if !haskey(kresults_dict, (x,z))
				kresults_dict[(x,z)] = []
			end
			push!(kresults_dict[(x,z)], kronresult)				
		end
	end

	result = Dict()
	for (label_pair, list_of_results) in kresults_dict
		
		buffer = Dict()		
		for kronresult in list_of_results
			for (newstates, value) in kronresult
				if !haskey(buffer, newstates)
					buffer[newstates]=S(0)
				end
				buffer[newstates]+=value
			end
		end
		result[label_pair] = buffer
	end
	result
	# cooC = SparseMatrixCOO{M, Int}(nsymbols, nsymbols, label_rows, label_cols, label_vals)
	# arcsC = coo_lod2arcs(cooC, Q, S)

	# initialA = zeros(S,numstates(A))
	# initialB = zeros(S,numstates(B))
	# initialA[A.initstate] = S(1)
	# initialB[B.initstate] = S(1)
	# initial = findfirst(x->x!=zero(S),kron(initialA, initialB))
	# final = kron(A.finalweights, B.finalweights) 

	# VectorFST(arcsC, initial, final )
	
end

# ╔═╡ 0199fbc6-44a7-48a6-a48c-0cee07aadfe8
begin
	cooA = dict2coo(vector2dict_lod(A), nsymbols, numstates(A), S)
	cooB = dict2coo(vector2dict_lod(B), nsymbols, numstates(B), S)
	Q = numstates(A) * numstates(B)
end

# ╔═╡ 07b0851b-c823-47b3-b6bf-921a0ed2bae2
@benchmark my2_sparse_composition_lod(cooA, cooB, S, Q, nsymbols)

# ╔═╡ df15d9c4-371d-45cf-8123-bb6f2d798932


# ╔═╡ 94869ac7-0935-47dd-b328-4decd3da01b1
S = TropicalSemiring{Float32}

# ╔═╡ 03b8a8cf-7677-44f9-a033-b991605ab297
# ╠═╡ disabled = true
#=╠═╡
S = ProbSemiring{Float32}
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═fdf50f2e-1ef7-11ee-0d03-db187df6536b
# ╠═03b8a8cf-7677-44f9-a033-b991605ab297
# ╠═b5f54187-2e6e-49f9-b327-6bfe501708c7
# ╠═dca59cae-aa9e-4dce-92b5-c4a2cca199ac
# ╠═f464e7ab-a831-473f-ad07-7f2de644956c
# ╠═8fadec26-cb9e-4533-ad8f-fffb376fd170
# ╠═46a6da30-8b40-4b94-a7a1-e01a428813cf
# ╠═94869ac7-0935-47dd-b328-4decd3da01b1
# ╠═a239c43c-9c10-406f-8f7b-c1b7d383dd3a
# ╠═145edb8f-9ffa-4e95-ba48-6c7dcdfeecd1
# ╠═efc53914-1c4e-4195-ac63-553f020c5f9b
# ╠═db03e5e8-5e31-4f7d-88f6-e3096f56b5c4
# ╠═2513f1a1-92cd-44cb-afe9-a6ef5dddffb1
# ╠═5c2e4fda-4304-4e6d-b98a-b51053f4b8ef
# ╠═4b777a57-8437-4486-b6e4-184c6801f1ae
# ╠═3e4c65ed-8ab8-4d88-8fef-24270f9946a7
# ╠═0199fbc6-44a7-48a6-a48c-0cee07aadfe8
# ╠═07b0851b-c823-47b3-b6bf-921a0ed2bae2
# ╠═df15d9c4-371d-45cf-8123-bb6f2d798932
