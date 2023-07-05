### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ b557cfde-1594-11ee-0a4a-83034cdfb96f
begin
	using Pkg
	Pkg.add(	
		[
			PackageSpec(name="Revise"),			
			PackageSpec(name="JSON"),
		]
	)

	Pkg.develop(path="finitestateautomata.jl/")

	using Revise
	using FiniteStateAutomata	
	using JSON
	using LinearAlgebra
	#using SparseSemimodules
end

# ╔═╡ 0ac086dc-c8a8-4584-bc3b-f987ae496db3


# ╔═╡ a33ebcad-67ed-4d95-b6b2-1443358f4e22
function diagm_(v)
	FiniteStateAutomata.sparse(collect(1:length(v)), collect(1:length(v)), v)
end

# ╔═╡ e906a5bd-2c8f-46ca-9830-659924d06773
STORAGE = "/Users/laneskij/Documents/phd-work/research/JSALT-2023/FAST/lecture/FST_JSALT23/libri_examples"

# ╔═╡ 97be68b2-1fd4-4a82-990b-25215f5e7c3e
begin
	label_mapping = open("$(STORAGE)/ctc_map.json") do f
		mapping = JSON.parse(f)
		Dict(parse(Int, k) + 1 => v for (k, v) in mapping)
	end
	label_mapping[length(label_mapping) + 1] = "BOS"
	label_mapping[length(label_mapping) + 1] = "EOS"
	label_mapping
end

# ╔═╡ d43b6234-a81c-42f2-98d6-8081ffa385cb
S = ProbSemiring{Float32}
#S = TropicalSemiring{Float32}

# ╔═╡ 4f4eadba-e495-48ec-9764-4565fb0a6894
struct DenseFST{S,L} <: FiniteStateAutomata.AbstractFST{S,L}
    M::Vector{Matrix{S}}
    α::Vector{S}
    ω::Vector{S}
    λ::AbstractVector{L}
end

# ╔═╡ e06ae027-6931-4a3d-a261-70cd42172fc2
begin
FiniteStateAutomata.M(fst::DenseFST) = fst.M
FiniteStateAutomata.α(fst::DenseFST) = fst.α
FiniteStateAutomata.ω(fst::DenseFST) = fst.ω
FiniteStateAutomata.λ(fst::DenseFST) = fst.λ	
end

# ╔═╡ 56b6ef39-c10f-4ce9-a87a-7462b29c97fa

A = DenseFST(
	[
		Array(sparse([1], [1], S[3], 4, 4)),
		Array(sparse([1, 2], [2, 3], S[1, 2], 4, 4)),
		Array(sparse([1, 2], [2, 3], S[3, 4], 4, 4)),		
		Array(sparse([1], [2], S[1], 4, 4)),
		Array(sparse([1, 3], [2, 3], S[5, 6], 4, 4)),
		Array(sparse([1, 3], [2, 4], S[7, 8], 4, 4)),
		Array(sparse([1, 2], [3, 4], S[9, 10], 4, 4)),
		
		
	],
	Array(sparsevec([1], one(S), 4)),
	Array(sparsevec([4], one(S), 4)),
	[0, 1, 6, 8, 13, 24, 2]
)

# ╔═╡ 9acc17d1-46d8-4b8d-b687-20ca90a3e303
draw(
	A;
	symbols = label_mapping
);

# ╔═╡ 0dd999dd-52f6-4696-9468-d6fff0b0370e
function NormFact(A::FiniteStateAutomata.AbstractFST)
	S = semiring(A)
	acc = deepcopy(A.ω)
	for n = 1:nstates(A)
		for m = 1:size(A.M)[1]
			acc[n] = acc[n] ⊕ sum(A.M[m][n,:])
		end
		acc[n] = inv(acc[n])
	end
		
acc
end
	
		

# ╔═╡ d1a62c1c-fdc5-4b83-bc9f-9e0d22d123d0
function LocallyNormalize(A::FiniteStateAutomata.AbstractFST, norm_fact::Array{S})
	ms = Matrix{S}[]
	for m = 1:size(A.M)[1]
		tmp = diagm(norm_fact) * A.M[m]
		@show typeof(tmp)
		push!(ms, tmp)
	end
	DenseFST(
		ms,
		A.α,
		A.ω,
		A.λ
	)
end

# ╔═╡ 521fcbf1-8aab-4282-8c6f-70bd9a2429f6


# ╔═╡ 8569b35b-886d-4b6e-9064-4e383aff2264
function get_transposed_matrix(A::FiniteStateAutomata.AbstractFST)
	S = semiring(A)
	
	ms = Matrix{S}[]
	for m = 1:size(A.M)[1]
		tmp = A.M[m].T
		push!(ms, tmp)
	end
	DenseFST(
		ms,
		A.α,
		A.ω,
		A.λ
	)
	
end

# ╔═╡ 0347b5ca-2bb4-4e2c-b5e0-970f61c1772c
trans_A = get_transposed_matrix(A)

# ╔═╡ aff61773-d2e1-4055-958e-c7b6299ab00b


# ╔═╡ c334ab7c-d31d-4351-b93b-291bd08f183b
# ╠═╡ disabled = true
#=╠═╡
nor
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═0ac086dc-c8a8-4584-bc3b-f987ae496db3
# ╠═a33ebcad-67ed-4d95-b6b2-1443358f4e22
# ╠═b557cfde-1594-11ee-0a4a-83034cdfb96f
# ╠═e906a5bd-2c8f-46ca-9830-659924d06773
# ╠═97be68b2-1fd4-4a82-990b-25215f5e7c3e
# ╠═d43b6234-a81c-42f2-98d6-8081ffa385cb
# ╠═4f4eadba-e495-48ec-9764-4565fb0a6894
# ╠═e06ae027-6931-4a3d-a261-70cd42172fc2
# ╠═56b6ef39-c10f-4ce9-a87a-7462b29c97fa
# ╠═9acc17d1-46d8-4b8d-b687-20ca90a3e303
# ╠═0dd999dd-52f6-4696-9468-d6fff0b0370e
# ╠═d1a62c1c-fdc5-4b83-bc9f-9e0d22d123d0
# ╠═521fcbf1-8aab-4282-8c6f-70bd9a2429f6
# ╠═8569b35b-886d-4b6e-9064-4e383aff2264
# ╠═0347b5ca-2bb4-4e2c-b5e0-970f61c1772c
# ╠═aff61773-d2e1-4055-958e-c7b6299ab00b
# ╠═c334ab7c-d31d-4351-b93b-291bd08f183b
