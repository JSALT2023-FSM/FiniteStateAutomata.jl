### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 6bf1296a-cfa0-11ed-0bfd-35ba2b7ac219
begin
	using Pkg
	Pkg.activate("..")

	using Revise
	using FiniteStateAutomata
end

# ╔═╡ f4e7e66b-71db-4a2a-a204-41c65a12062b
K = Float32

# ╔═╡ 3a3deafb-4154-4abf-b7b6-3c0317cd7457
H = randn(K, 3, 4)

# ╔═╡ f85a5dd1-9e82-40cf-a843-df2a658bc867
DenseFSA(H, ["a", "b", "c"])

# ╔═╡ Cell order:
# ╠═6bf1296a-cfa0-11ed-0bfd-35ba2b7ac219
# ╠═f4e7e66b-71db-4a2a-a204-41c65a12062b
# ╠═3a3deafb-4154-4abf-b7b6-3c0317cd7457
# ╠═f85a5dd1-9e82-40cf-a843-df2a658bc867
