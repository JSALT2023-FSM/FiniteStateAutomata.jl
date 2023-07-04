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

# ╔═╡ b64e8b2b-9c69-4cce-9595-a04265c9a299
md"
## FST Composition Example
### States first ordering
"

# ╔═╡ 1f665bf9-99a2-4244-a0a8-87d5c8a1fcb0
begin
	C = sparse_composition_kron(A, B, length(symbols))
	draw(C; symbols)
end

# ╔═╡ Cell order:
# ╟─38eca186-37ed-4ae7-8341-3e28fee6ecc0
# ╠═7911b2c4-183f-11ee-05e9-c9b311ee4780
# ╠═2c3c6583-491c-4af9-a424-7d5f983ee982
# ╠═438a3823-a65f-4ab3-9cf0-7ebf946c9226
# ╠═fd4f5379-7d6b-4962-8f07-ddedacf565a6
# ╠═ddff62e9-b7ed-43b3-beaf-c46c8ebe61ec
# ╠═2d6277aa-d73f-43be-9a91-2e82218a5587
# ╟─b64e8b2b-9c69-4cce-9595-a04265c9a299
# ╠═1f665bf9-99a2-4244-a0a8-87d5c8a1fcb0
