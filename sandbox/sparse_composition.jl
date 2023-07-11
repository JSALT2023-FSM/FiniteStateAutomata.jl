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
draw(A, isymbols=symbols, osymbols=symbols) |> Dot2SVG() |> HTML

# ╔═╡ 2d6277aa-d73f-43be-9a91-2e82218a5587
draw(B, isymbols=symbols, osymbols=symbols) |> Dot2SVG() |> HTML

# ╔═╡ c12bb3ec-c6cb-4467-b221-84cb812ff549
draw(B, isymbols=symbols, osymbols=symbols) |> Dot2SVG() 

# ╔═╡ b64e8b2b-9c69-4cce-9595-a04265c9a299
md"
## FST Composition Example
### States first ordering
"

# ╔═╡ 1f665bf9-99a2-4244-a0a8-87d5c8a1fcb0
begin
	C = sparse_composition_sod(A, B, length(symbols))
	draw(C, isymbols=symbols, osymbols=symbols) |> Dot2SVG() |> HTML
end

# ╔═╡ ee5d158a-5855-4b4d-be2f-591139ee487d
Vector{Tuple{Int64, Int64, Int64, FiniteStateAutomata.TropicalSemiring{Float64}}}[[(5, 1, 3, 0.4)], [(9, 2, 2, 0.6000000000000001)], [(9, 2, 2, 0.8)], [], [(12, 1, 2, 0.8), (6, 3, 2, 0.7)], [(12, 1, 2, 1.0), (6, 3, 2, 0.8999999999999999)], [(11, 2, 3, 0.8)], [], [], [], [], []]

# ╔═╡ c5c8ecf4-db34-4146-bf60-22906abb452a
md"
### Labels first ordering
"

# ╔═╡ d762fb36-1c9c-4f1d-9e2d-683c74b92693
begin   
	C2 = sparse_composition_lod(A, B, length(symbols))
	draw(C2, isymbols=symbols, osymbols=symbols) |> Dot2SVG() |> HTML
end 

# ╔═╡ 6ef546a1-a2d1-4a92-a190-6fd6eb43a59c
C2.arcs==C.arcs

# ╔═╡ Cell order:
# ╟─38eca186-37ed-4ae7-8341-3e28fee6ecc0
# ╠═7911b2c4-183f-11ee-05e9-c9b311ee4780
# ╠═2c3c6583-491c-4af9-a424-7d5f983ee982
# ╠═438a3823-a65f-4ab3-9cf0-7ebf946c9226
# ╠═fd4f5379-7d6b-4962-8f07-ddedacf565a6
# ╠═ddff62e9-b7ed-43b3-beaf-c46c8ebe61ec
# ╠═2d6277aa-d73f-43be-9a91-2e82218a5587
# ╠═c12bb3ec-c6cb-4467-b221-84cb812ff549
# ╟─b64e8b2b-9c69-4cce-9595-a04265c9a299
# ╠═1f665bf9-99a2-4244-a0a8-87d5c8a1fcb0
# ╠═ee5d158a-5855-4b4d-be2f-591139ee487d
# ╟─c5c8ecf4-db34-4146-bf60-22906abb452a
# ╠═d762fb36-1c9c-4f1d-9e2d-683c74b92693
# ╠═6ef546a1-a2d1-4a92-a190-6fd6eb43a59c
