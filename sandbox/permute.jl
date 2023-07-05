### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ a183e6f2-157e-11ee-0d03-ad978a190e29
begin
	using Pkg
	Pkg.activate("../")

	using Images
	using Revise
	using FiniteStateAutomata
end

# ╔═╡ dcbb87c3-139b-465c-a1e0-04774a868411
S = LogSemiring{Float32,1}

# ╔═╡ da6fdaa9-675d-49b8-80e7-a5aaae988b10
symsdir = "../examples/symboltables/"

# ╔═╡ 49fead47-9d06-4919-a675-ea1591aa83eb
syms = open(loadsymbols, "../examples/symboltables/ascii.syms")

# ╔═╡ 96d478e3-f21f-4fa8-ad76-3f88424d3a52
fst = open("../examples/fsts/fst_ex1.txt") do f
		compile(f; semiring = S)
	end

# ╔═╡ 5a738e4d-5667-468a-8090-4f94ca4a413c
display(MIME("image/svg+xml"), fst)

# ╔═╡ 70d3f948-70b7-4a40-8346-0d190e70f77d
x = draw(fst)

# ╔═╡ 57395d0d-567d-4329-9794-2a2a1b525efa
display(MIME("image/png"), draw(fst))

# ╔═╡ 0cdf29b1-09db-43ca-8aab-b7e60f4928ad
HTML("""
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
 "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<!-- Generated by graphviz version 8.0.5 (20230430.1635)
 -->
<!-- Pages: 1 -->
<svg width="300pt" height="78pt"
 viewBox="0.00 0.00 299.91 77.87" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
<g id="graph0" class="graph" transform="scale(1 1) rotate(0) translate(4 73.87)">
<polygon fill="white" stroke="none" points="-4,4 -4,-73.87 295.91,-73.87 295.91,4 -4,4"/>
<!-- 1 -->
<g id="node1" class="node">
<title>1</title>
<ellipse fill="none" stroke="black" stroke-width="2" cx="18.38" cy="-34.94" rx="18.38" ry="18.38"/>
<text text-anchor="middle" x="18.38" y="-30.64" font-family="Times,serif" font-size="14.00">1</text>
</g>
<!-- 2 -->
<g id="node2" class="node">
<title>2</title>
<ellipse fill="none" stroke="black" cx="129.4" cy="-34.94" rx="18.38" ry="18.38"/>
<text text-anchor="middle" x="129.4" y="-30.64" font-family="Times,serif" font-size="14.00">2</text>
</g>
<!-- 1&#45;&gt;2 -->
<g id="edge1" class="edge">
<title>1&#45;&gt;2</title>
<path fill="none" stroke="black" d="M36.9,-34.94C53.81,-34.94 79.7,-34.94 99.63,-34.94"/>
<polygon fill="black" stroke="black" points="99.57,-38.44 109.57,-34.94 99.57,-31.44 99.57,-38.44"/>
<text text-anchor="middle" x="73.89" y="-39.64" font-family="Times,serif" font-size="14.00">1:1/0.5</text>
</g>
<!-- 1&#45;&gt;2 -->
<g id="edge2" class="edge">
<title>1&#45;&gt;2</title>
<path fill="none" stroke="black" d="M33.32,-23.76C39.51,-19.49 47.12,-15.12 54.77,-12.94 71.12,-8.27 76.67,-8.27 93.02,-12.94 97.2,-14.13 101.38,-15.98 105.31,-18.11"/>
<polygon fill="black" stroke="black" points="103.12,-21.49 113.47,-23.76 106.8,-15.53 103.12,-21.49"/>
<text text-anchor="middle" x="73.89" y="-17.64" font-family="Times,serif" font-size="14.00">2:2/1.5</text>
</g>
<!-- 3 -->
<g id="node3" class="node">
<title>3</title>
<ellipse fill="none" stroke="black" cx="256.98" cy="-34.94" rx="30.94" ry="30.94"/>
<ellipse fill="none" stroke="black" cx="256.98" cy="-34.94" rx="34.94" ry="34.94"/>
<text text-anchor="middle" x="256.98" y="-30.64" font-family="Times,serif" font-size="14.00">3/3.5</text>
</g>
<!-- 2&#45;&gt;3 -->
<g id="edge3" class="edge">
<title>2&#45;&gt;3</title>
<path fill="none" stroke="black" d="M148.13,-34.94C164.44,-34.94 189.32,-34.94 211.08,-34.94"/>
<polygon fill="black" stroke="black" points="210.78,-38.44 220.78,-34.94 210.78,-31.44 210.78,-38.44"/>
<text text-anchor="middle" x="184.91" y="-39.64" font-family="Times,serif" font-size="14.00">3:3/1.5</text>
</g>
</g>
</svg>
""")

# ╔═╡ Cell order:
# ╠═a183e6f2-157e-11ee-0d03-ad978a190e29
# ╠═dcbb87c3-139b-465c-a1e0-04774a868411
# ╠═da6fdaa9-675d-49b8-80e7-a5aaae988b10
# ╠═49fead47-9d06-4919-a675-ea1591aa83eb
# ╠═96d478e3-f21f-4fa8-ad76-3f88424d3a52
# ╠═5a738e4d-5667-468a-8090-4f94ca4a413c
# ╠═70d3f948-70b7-4a40-8346-0d190e70f77d
# ╠═57395d0d-567d-4329-9794-2a2a1b525efa
# ╠═0cdf29b1-09db-43ca-8aab-b7e60f4928ad