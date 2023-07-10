### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 1776376e-1bfd-11ee-07bc-dfd83b6e91f0
begin
	using Pkg
	Pkg.develop(path="../../finitestateautomata.jl/")	
	using Revise
	using FiniteStateAutomata
    using PlutoUI, BenchmarkTools
	using Profile
	using ProfileCanvas
	using NPZ
end

# ╔═╡ 2c64ec0d-406b-4886-8e4d-8dc9caa22622
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
</div>""")

# ╔═╡ 6610af7b-f81b-421f-9bb3-49d44bfc9e60
S = TropicalSemiring{Float32}

# ╔═╡ 50ceeb69-1dbc-447f-8584-de72d9a08473
symbols = Dict(1 => "a", 2 => "b", 3 => "c")

# ╔═╡ 5e310b5f-49f6-485e-afa4-8b4ae66cd56b
A = convert(TensorFST{S, Array{S,4}}, VectorFST(
	[
		Arc{S}[(2, 1, 2, S(0.1)), (3, 2, 1, S(0.2))],
		Arc{S}[(2, 3, 1, S(0.3)), (4, 1, 1, S(0.4))],
		Arc{S}[(4, 2, 2, S(0.5))],
		Arc{S}[]
	],
	1,
	[zero(S),zero(S),zero(S),S(0.6)]
));

# ╔═╡ 88353a05-4f73-42c3-9245-11fa9a8360a8
B = convert(TensorFST{S, Array{S,4}}, VectorFST(
	[
		Arc{S}[(2, 2, 3, S(0.3))],
		Arc{S}[(3, 1, 2, S(0.4))],
		Arc{S}[(3, 1, 2, S(0.6))],
	],
	1,
	S[zero(S),zero(S),S(0.7)]
));

# ╔═╡ 9453ddc9-cdda-4aa6-ae1c-51392639b5a4
@benchmark dense_composition_sod(A, B)

# ╔═╡ 0da10c31-7c3a-4560-b5d5-6e05b43cf084
@benchmark dense_composition_lod(A, B)

# ╔═╡ 0e1e2f9d-fa73-440e-a4eb-ef639babbcd7
md"## Loading an ASR FST"

# ╔═╡ b6ffd137-5ff4-448a-8a6e-0a0747b19454
logits = NPZ.npzread("../examples/assets/libri_examples/2830-3980-0002/logits.npy");

# ╔═╡ 9985a493-cbfd-4b9c-86dc-8e1c080228bb
size(logits)

# ╔═╡ 20ddef92-0e66-4c9a-b166-ec3b91035c6a
begin
	label_mapping = open("../examples/assets/libri_examples/ctc_map.txt") do f
		readlines(f)
	end
	label_mapping = Dict(k => v for (k,v) in enumerate(label_mapping))
	label_mapping_reverse = Dict(value => key for (key, value) in label_mapping)
	nsym = length(label_mapping)
end;

# ╔═╡ b75b944e-cf79-4268-9891-07f95af4d43c
label_mapping

# ╔═╡ 24189d1c-83c6-4555-ab1c-7d15a2eab134
begin
	# lets represent the logits as label-dependent alignment lattice
	local Q, L = size(logits)
	Q += 1 # +1 for initial/final state
	local K = LogSemiring{Float16, ℯ}
	
	local W = K.(logits)
	local α = zeros(K, Q)
	α[1] = one(K)
	local ω = zeros(K, Q)
	ω[end] = one(K)
	local M = zeros(K, Q, Q, L, L)

	for q in 1:(Q-1)
		for l in 1:L
			M[q, q+1, l, l] = W[q, l]
 		end
	end
	nn_tensor_fst = TensorFST(M, α, ω)


	local arc_type = Vector{Tuple{Int,Int,Int,K}}		
	 
	local myarcs = Vector{arc_type}()
	for i in 1:Q
		push!(myarcs,Vector{arc_type}())
	end

	for q in 1:(Q-1)
		for l in 1:L
			push!(myarcs[q], (q+1, l, l, W[q, l] ))
 		end
	end
	
	local final = zeros(K,Q)
	final[Q] = K(1.0)
	nn_vector_fst = VectorFST(myarcs, 1, final)
	nn_sparse_sod_fst = dict2coo(vector2dict_sod(nn_vector_fst), Q, nsym, K)
	nn_sparse_lod_fst = dict2coo(vector2dict_lod(nn_vector_fst), Q, nsym, K)
end;

# ╔═╡ c657cae1-b442-4f65-a27b-5738c20dad8e
# Linear FST from word in vector format
"""
Word as sequence of character\n
char2index mapping\n
K semiring\n
kw weight for repeating state with blank\n
nw weight for next character\n
"""
function word2linearwloops_fst(word, char2index, K, kw, nw)
	Sprob = ProbSemiring{Float32}
	local arc_type = Vector{Tuple{Int,Int,Int,K}}		
	Q = length(word)+1
	myarcs = Vector{arc_type}()
	for i in 1:Q
		push!(myarcs,Vector{arc_type}())
	end
	for i in 1:Q-1
		k = char2index[string(word[i])]
		push!(myarcs[i], (i, 1 , 1, kw ))
		push!(myarcs[i], (i+1, k , k, nw ))
	end
	push!(myarcs[Q], (Q, 1 , 1, kw ))

	local final = zeros(K,Q)
	final[Q] = K(1.0)
	VectorFST(myarcs, 1, final)
end

# ╔═╡ 55a29b8b-2153-40f6-b5f5-47bf55bfc822
begin
	K = LogSemiring{Float16, ℯ}
	word_vector_fst = word2linearwloops_fst("HELLO", label_mapping_reverse, K, K(0.1), K(1))
	draw(word_vector_fst, isymbols=label_mapping, osymbols=label_mapping) |> Dot2SVG() |> HTML
end

# ╔═╡ 0dedbc99-bdb9-425f-a199-3d636184e741
begin
	comp_sod = sparse_composition_sod(nn_vector_fst, word_vector_fst, nsym)
	draw(comp_sod, isymbols=label_mapping, osymbols=label_mapping) |> Dot2SVG() |> HTML
end

# ╔═╡ 3096cc0a-0e69-459c-a735-f8a209b0adf6
begin
	comp_lod = sparse_composition_lod(nn_vector_fst, word_vector_fst, nsym)
	draw(comp_lod, isymbols=label_mapping, osymbols=label_mapping) |> Dot2SVG() |> HTML
end

# ╔═╡ 04d1e73f-fdb7-4e68-a5c9-438d4b7bb15f
@benchmark sparse_composition_lod(nn_vector_fst, word_vector_fst, nsym)

# ╔═╡ 09e889f8-1a0b-47ed-abcc-3907625ec417
@benchmark sparse_composition_sod(nn_vector_fst, word_vector_fst, nsym)

# ╔═╡ 32081a1f-f442-45b0-8c22-b2ea47a99ffb
"States",numstates(nn_vector_fst), "arcs", sum(map(length,nn_vector_fst.arcs))

# ╔═╡ 696a18fb-0619-4069-8699-e32a3be22819
begin
	charlm = open(x-> compile(x,semiring=K, openfst_compat=true), "../../notebooks/charwotw.2.fst.txt") 
	# setinitstate!(lm,1) 
end

# ╔═╡ a4a1e489-bb36-4589-b8d1-84578ee6b39d
"States",numstates(charlm), "arcs", sum(map(length,charlm.arcs))

# ╔═╡ f534b0c8-af0d-4fda-b2cf-9da726af36d3
# draw(lm, isymbols=label_mapping, osymbols=label_mapping) |> Dot2PNG() |> HTML

# ╔═╡ e930e517-66a0-4e77-8962-db9c0d11b2e5
# open("../../notebooks/ctc_map.txt", "r") do f
# 	print(read(f, String))
# end

# ╔═╡ 69e3408e-5f42-4a84-8174-29ed35b8ce7d
@benchmark sparse_composition_sod(nn_vector_fst, charlm, nsym)

# ╔═╡ f98b0a0a-7194-45ee-b43e-c1de9e59de38
@benchmark sparse_composition_lod(nn_vector_fst, charlm, nsym)

# ╔═╡ 5f1494f6-372a-414e-b960-5cd8d51d3b7e
@benchmark sparse_composition_lod_mt(nn_vector_fst, charlm, nsym)   

# ╔═╡ 5cb7529a-e9f5-4110-8730-63a8541017c3
md"Bench with out conversions"

# ╔═╡ b101405c-22a6-400b-8c33-3a2b3453f560
begin
	local A = nn_vector_fst
	local B = charlm
	local S = semiring(A)
	Q = numstates(A) * numstates(B)
	cooAsod = dict2coo(vector2dict_sod(A), numstates(A),  nsym, S)
	cooBsod = dict2coo(vector2dict_sod(B), numstates(B),  nsym, S)
end;

# ╔═╡ a078eaf3-22c3-489f-9f18-e837c7aeea8b
@benchmark cooC = kron(cooAsod, cooBsod)	

# ╔═╡ 58a9f075-0266-4d98-a25a-f3726905d019
begin
	local A = nn_vector_fst
	local B = charlm
	local S = semiring(A)
	cooAlod = dict2coo(vector2dict_lod(A), nsym, numstates(A), S)
	cooBlod = dict2coo(vector2dict_lod(B), nsym, numstates(B), S)
end;

# ╔═╡ 1306d047-5e73-498e-a37c-d6a9207653ad
# sparse_coo2dict_composition_lod_mt(cooAlod, cooBlod, semiring(nn_vector_fst), Q, nsym)

# ╔═╡ 965a8d1a-d909-4f83-a334-205e703ffcbf
@benchmark sparse_coo2dict_composition_lod_mt(cooAlod, cooBlod, semiring(nn_vector_fst), Q, nsym)

# ╔═╡ 88e0f162-57a9-49ea-8b28-8b18430b4b4f
# begin
# 	function profile_test_sod(n)
# 		for i = 1:n
# 			 sparse_composition_sod(nn_vector_fst, charlm, nsym)
# 		end
# 	end
	
# 	function profile_test_lod(n)
# 		for i = 1:n
# 			 sparse_composition_lod(nn_vector_fst, charlm, nsym)
# 		end
# 	end
# end

# ╔═╡ 38b960c2-0bef-4069-bf5d-0b44649e4a25
# @profview profile_test_sod(1);

# ╔═╡ bca99001-10bb-44e8-ad7d-2f5a3f602502
# @profview profile_test_sod(10)

# ╔═╡ e62a648d-2667-4105-99aa-a8ba375f37d1
# @profview profile_test_lod(1);

# ╔═╡ 54253fbb-01bd-400c-bffb-5f62430093d7
# @profview profile_test_lod(10)

# ╔═╡ Cell order:
# ╟─2c64ec0d-406b-4886-8e4d-8dc9caa22622
# ╠═1776376e-1bfd-11ee-07bc-dfd83b6e91f0
# ╠═6610af7b-f81b-421f-9bb3-49d44bfc9e60
# ╠═50ceeb69-1dbc-447f-8584-de72d9a08473
# ╠═5e310b5f-49f6-485e-afa4-8b4ae66cd56b
# ╠═88353a05-4f73-42c3-9245-11fa9a8360a8
# ╠═9453ddc9-cdda-4aa6-ae1c-51392639b5a4
# ╠═0da10c31-7c3a-4560-b5d5-6e05b43cf084
# ╠═0e1e2f9d-fa73-440e-a4eb-ef639babbcd7
# ╠═b6ffd137-5ff4-448a-8a6e-0a0747b19454
# ╠═9985a493-cbfd-4b9c-86dc-8e1c080228bb
# ╠═20ddef92-0e66-4c9a-b166-ec3b91035c6a
# ╠═b75b944e-cf79-4268-9891-07f95af4d43c
# ╠═24189d1c-83c6-4555-ab1c-7d15a2eab134
# ╠═c657cae1-b442-4f65-a27b-5738c20dad8e
# ╠═55a29b8b-2153-40f6-b5f5-47bf55bfc822
# ╠═0dedbc99-bdb9-425f-a199-3d636184e741
# ╠═3096cc0a-0e69-459c-a735-f8a209b0adf6
# ╠═04d1e73f-fdb7-4e68-a5c9-438d4b7bb15f
# ╠═09e889f8-1a0b-47ed-abcc-3907625ec417
# ╠═32081a1f-f442-45b0-8c22-b2ea47a99ffb
# ╠═696a18fb-0619-4069-8699-e32a3be22819
# ╠═a4a1e489-bb36-4589-b8d1-84578ee6b39d
# ╠═f534b0c8-af0d-4fda-b2cf-9da726af36d3
# ╠═e930e517-66a0-4e77-8962-db9c0d11b2e5
# ╠═69e3408e-5f42-4a84-8174-29ed35b8ce7d
# ╠═f98b0a0a-7194-45ee-b43e-c1de9e59de38
# ╠═5f1494f6-372a-414e-b960-5cd8d51d3b7e
# ╟─5cb7529a-e9f5-4110-8730-63a8541017c3
# ╠═b101405c-22a6-400b-8c33-3a2b3453f560
# ╠═a078eaf3-22c3-489f-9f18-e837c7aeea8b
# ╠═58a9f075-0266-4d98-a25a-f3726905d019
# ╠═1306d047-5e73-498e-a37c-d6a9207653ad
# ╠═965a8d1a-d909-4f83-a334-205e703ffcbf
# ╠═88e0f162-57a9-49ea-8b28-8b18430b4b4f
# ╠═38b960c2-0bef-4069-bf5d-0b44649e4a25
# ╠═bca99001-10bb-44e8-ad7d-2f5a3f602502
# ╠═e62a648d-2667-4105-99aa-a8ba375f37d1
# ╠═54253fbb-01bd-400c-bffb-5f62430093d7
