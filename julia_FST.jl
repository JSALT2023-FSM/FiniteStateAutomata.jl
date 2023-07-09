### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ dc65255e-7f62-47ff-9b13-b303fbf6793c
begin 
	using Pkg
	Pkg.develop(path="julia_FST/finitestateautomata.jl/")
	using Zygote
	using FiniteStateAutomata
	using Zygote: @adjoint
end


# ╔═╡ a5ca177f-5116-408a-be1e-0469de45b1c6
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

# ╔═╡ dde7c4da-55f3-4e89-afd3-e81ae5fd33a7
typeof(LogSemiring{Float32, 2}(2))

# ╔═╡ 55dbf49d-564f-419c-aef0-1cf1951e38b9
LogSemiring2 = LogSemiring{Float32, 2}

# ╔═╡ 7c6fbecb-ebb5-4031-b0a6-09063f8e11d0
ProbSemiring(3.0)

# ╔═╡ 5dbcc108-724b-4f9c-92ee-929454eab873
g(x) = ProbSemiring(2.0) ⊗ ProbSemiring(x) ⊗ ProbSemiring(x)

# ╔═╡ cd821576-5cf1-4a32-9297-c5fb106ed5bf
q(x) = LogSemiring2(2.0) ⊗ LogSemiring2(x) ⊗ LogSemiring2(x)

# ╔═╡ c0a6c01c-17c9-4ab0-9d7b-d0e7cc388a9d
@adjoint ProbSemiring(x) = ProbSemiring(x), a -> (a,)

# ╔═╡ f90de6a6-2178-4b72-9b80-e12abd75e4a9
@adjoint LogSemiring{Float32, w}(x) = LogSemiring{Float32, w}(x), a -> (a,)

# ╔═╡ 603e34fc-f175-4626-b202-15704953f94c
@adjoint T(x) = T(x), a -> (a,) with T<:Semiring

# ╔═╡ d1165269-fe0c-42c2-b998-f85d479aadf6
q'(2.0)

# ╔═╡ f3fc8337-7a9e-45fd-87f9-803394fde566
g'(2.0) 

# ╔═╡ 4ef031b9-6bd9-4c72-a8f2-acd58baf3887
f(x) = x*x*5 + 2*x + 3

# ╔═╡ 8c02f5df-fea7-447d-923c-eae7964190a6
f'(10)

# ╔═╡ 4a4292d9-9564-450f-bc9e-652a457c99a7


# ╔═╡ 1cb35cd0-d7d1-4643-a7c4-d7ddf9a7e782


# ╔═╡ aeb11514-119c-4aa7-a0e8-22d6db86ecd7


# ╔═╡ Cell order:
# ╟─a5ca177f-5116-408a-be1e-0469de45b1c6
# ╠═dc65255e-7f62-47ff-9b13-b303fbf6793c
# ╠═dde7c4da-55f3-4e89-afd3-e81ae5fd33a7
# ╠═55dbf49d-564f-419c-aef0-1cf1951e38b9
# ╠═7c6fbecb-ebb5-4031-b0a6-09063f8e11d0
# ╠═5dbcc108-724b-4f9c-92ee-929454eab873
# ╠═cd821576-5cf1-4a32-9297-c5fb106ed5bf
# ╠═c0a6c01c-17c9-4ab0-9d7b-d0e7cc388a9d
# ╠═f90de6a6-2178-4b72-9b80-e12abd75e4a9
# ╠═603e34fc-f175-4626-b202-15704953f94c
# ╠═d1165269-fe0c-42c2-b998-f85d479aadf6
# ╠═f3fc8337-7a9e-45fd-87f9-803394fde566
# ╠═4ef031b9-6bd9-4c72-a8f2-acd58baf3887
# ╠═8c02f5df-fea7-447d-923c-eae7964190a6
# ╠═4a4292d9-9564-450f-bc9e-652a457c99a7
# ╠═1cb35cd0-d7d1-4643-a7c4-d7ddf9a7e782
# ╠═aeb11514-119c-4aa7-a0e8-22d6db86ecd7
