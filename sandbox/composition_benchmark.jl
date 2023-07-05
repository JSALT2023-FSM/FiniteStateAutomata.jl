### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ e01251a6-1a74-11ee-2d45-45597e2308c3
begin
	using Pkg
	Pkg.develop(path="../../finitestateautomata.jl/")	
	using Revise
	using FiniteStateAutomata
    using PlutoUI, BenchmarkTools
	using Profile
end

# ╔═╡ 46a194ec-2fe3-4f7d-8d1c-5fe7bd63a760
using NPZ

# ╔═╡ 7880a62f-d158-472f-8056-11e56bf0062f
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

# ╔═╡ 857abf0c-7eee-4d17-86da-4e647f4310f2
S = TropicalSemiring{Float32}

# ╔═╡ 9576d2db-d859-49d1-8f1d-2b9394fb0152
symbols = Dict(1 => "a", 2 => "b", 3 => "c")

# ╔═╡ aae816b6-c947-4c23-9a6b-397cee58af5b
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

# ╔═╡ 393b3311-8efa-404b-9c17-406f7208db9e
B = convert(TensorFST{S, Array{S,4}}, VectorFST(
	[
		Arc{S}[(2, 2, 3, S(0.3))],
		Arc{S}[(3, 1, 2, S(0.4))],
		Arc{S}[(3, 1, 2, S(0.6))],
	],
	1,
	S[zero(S),zero(S),S(0.7)]
));

# ╔═╡ 3b3cb9f2-84e7-472b-836e-ec397a93166c
@benchmark dense_composition_sfo(A, B)

# ╔═╡ 18a6a191-605f-4705-b1eb-60b6421c5b9c
@benchmark dense_composition_lfo(A, B)

# ╔═╡ a730886c-ab35-4008-8069-452cab13f9cf
md"## Loading an ASR FST"

# ╔═╡ 1fe61a92-c18a-4cad-acdb-acf30dd36b7d
logits = NPZ.npzread("../examples/assets/libri_examples/2830-3980-0002/logits.npy");

# ╔═╡ 7eaae0a1-0300-4dcd-b3d6-104f6d182b86
begin
	label_mapping = open("../examples/assets/libri_examples/ctc_map.txt") do f
		readlines(f)
	end
	label_mapping = Dict(k => v for (k,v) in enumerate(label_mapping))
	label_mapping_reverse = Dict(value => key for (key, value) in label_mapping)
	nsym = length(label_mapping)
end;

# ╔═╡ 19692dee-38d9-4584-8adc-6b3075608269
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
	nn_sparse_fst = dict2coo(vector2dict(nn_vector_fst), Q, nsym, K)		
end;

# ╔═╡ 5c086ce1-1872-4750-b3cd-6e655a0e1a2e
begin
	local K = LogSemiring{Float16, ℯ}
	word = "HELLO"
	kw = K(0.1)
	nw = K(1)
	Sprob = ProbSemiring{Float32}
	local arc_type = Vector{Tuple{Int,Int,Int,K}}		
	Q = length(word)+1
	myarcs = Vector{arc_type}()
	for i in 1:Q
		push!(myarcs,Vector{arc_type}())
	end
	for i in 1:Q-1
		k = label_mapping_reverse[string(word[i])]
		push!(myarcs[i], (i, 1 , 1, kw ))
		push!(myarcs[i], (i+1, k , k, nw ))
	end
	push!(myarcs[Q], (Q, 1 , 1, kw ))

	local final = zeros(K,Q)
	final[Q] = K(1.0)
	word_vector_fst = VectorFST(myarcs, 1, final)
	word_sparse_fst = dict2coo(vector2dict(word_vector_fst), Q, nsym, K)		

	T = zeros(K, Q, Q, nsym, nsym) # shape Q x Q x L x L
    for (src, arcs) in enumerate(word_vector_fst.arcs)
        for (dest, isym, osym, w) in arcs
            T[src, dest, isym, osym] = w
        end
    end
    α = zeros(K, numstates(word_vector_fst))
    α[word_vector_fst.initstate] = one(K)
    ω = word_vector_fst.finalweights
    word_tensor_fst = TensorFST(T, α, ω)
	"Creates Linear FST from word in TensorFST format"
end

# ╔═╡ 56cbef6f-f4fe-49d2-b87f-f1b54569164a
draw(word_vector_fst, symbols = label_mapping)

# ╔═╡ b87c3d9d-be46-47c1-b215-c28ffbe8ece0
# nn_word_fst = dense_composition_lfo(nn_tensor_fst, word_tensor_fst);
# nn_word_fst2 = dense_composition_sfo(nn_fst, word_fst);

# ╔═╡ 3820427c-2c5b-49bc-b909-56ba0858c72b
# draw(nn_word_fst, symbols=label_mapping)

# ╔═╡ de4afc37-2f02-4830-a9f1-6781bb958292
comp = sparse_composition_sfo(nn_vector_fst, word_vector_fst, nsym)

# ╔═╡ 1c144d3f-85e8-4d35-b34e-a24da9e63b0e
draw(comp, symbols=label_mapping)

# ╔═╡ Cell order:
# ╟─7880a62f-d158-472f-8056-11e56bf0062f
# ╠═e01251a6-1a74-11ee-2d45-45597e2308c3
# ╠═857abf0c-7eee-4d17-86da-4e647f4310f2
# ╠═9576d2db-d859-49d1-8f1d-2b9394fb0152
# ╠═aae816b6-c947-4c23-9a6b-397cee58af5b
# ╠═393b3311-8efa-404b-9c17-406f7208db9e
# ╠═3b3cb9f2-84e7-472b-836e-ec397a93166c
# ╠═18a6a191-605f-4705-b1eb-60b6421c5b9c
# ╟─a730886c-ab35-4008-8069-452cab13f9cf
# ╠═46a194ec-2fe3-4f7d-8d1c-5fe7bd63a760
# ╠═1fe61a92-c18a-4cad-acdb-acf30dd36b7d
# ╠═7eaae0a1-0300-4dcd-b3d6-104f6d182b86
# ╠═19692dee-38d9-4584-8adc-6b3075608269
# ╠═5c086ce1-1872-4750-b3cd-6e655a0e1a2e
# ╠═56cbef6f-f4fe-49d2-b87f-f1b54569164a
# ╠═b87c3d9d-be46-47c1-b215-c28ffbe8ece0
# ╠═3820427c-2c5b-49bc-b909-56ba0858c72b
# ╠═de4afc37-2f02-4830-a9f1-6781bb958292
# ╠═1c144d3f-85e8-4d35-b34e-a24da9e63b0e
