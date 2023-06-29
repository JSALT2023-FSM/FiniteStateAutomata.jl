### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 0a906f46-1ee0-4814-860a-56e127b71593
begin
	using Pkg
	Pkg.develop(path="../repos/finitestateautomata.jl/")
	using Revise
	using FiniteStateAutomata
    using PlutoUI, BenchmarkTools
	using ProfileCanvas
	using Profile
end

# ╔═╡ 8b30c1e2-f00d-4d8c-b91a-9d7401a4628f
S = TropicalSemiring{Float32}

# ╔═╡ 7e5488a5-4de3-4543-91f5-09e8db4f5c36
symbols = Dict(1 => "a", 2 => "b", 3 => "c")

# ╔═╡ 314ba47a-2159-42ab-b580-f55759c75169
function external_dot(S, x,y)
	# dummy external product
	z = zeros(S,(length(x), length(y)))
	for (i,vx) in enumerate(x)
		for (j,vy) in enumerate(y)
			if vx!=zero(S) && vy!=zero(S)
				@show i,j, vx, vy
				z[i,j] = vx*vy
			end
		end
	end
	z
end

# ╔═╡ 5723bafe-c3e1-4b8f-8d22-7f9d6292f7f1
function coo2lil(A)
	t = Vector{Tuple{Int,Pair{Int,Int},S}}
	nstates = maximum(maximum(A,dims=1)[1][1:2])
	states = Vector{t}()
	for i in 1:nstates
		push!(states,Vector{t}())
	end
	for (i,j,k,l,v) in A
		push!(states[i],(j,k=>l,v))
	end	
	states
end

# ╔═╡ 699a88dc-ec13-4a4a-bdc9-f554a4be5e61
function dense2coo(A)
	coords = findall(!iszero,A)
	arcs = []
	for c in coords
		i,j,k,l = Tuple(c)
	 	push!(arcs,[i,j,k,l,A[i,j,k,l]])
	end
	arcs
end

# ╔═╡ 807cde31-89ca-4b45-984d-f80d01ace895
function coo2dense(S,A,nstates,nsymbols)
	TA = zeros(S,(nstates,nstates,nsymbols,nsymbols))
	for (i,j,k,l,v) in A
		TA[i,j,k,l] = v
	end
	TA
end

# ╔═╡ 08e20858-4f43-4614-9d64-fcf7a32c7fa8
function composition(DA,DB, nsymbols)
	n = size(DA,1)
	m = size(DB,1)
	DT = zeros(S,(n*m,n*m,nsymbols,nsymbols))
	for i in 1:n
		for j in 1:n
			for k in 1:m
				for l in 1:m
					DT[(k-1)*n+i,(l-1)*n+j,:,:] .= DA[i,j,:,:]*DB[k,l,:,:]
				end
			end
		end
	end
	DT
end

# ╔═╡ c7b97138-1674-11ee-1f7b-e9882a147ade
begin
	A = [
		[1,2,1,2,S(0.1)],
		[1,3,2,1,S(0.2)],
		[2,2,3,1,S(0.3)],
		[2,4,1,1,S(0.4)],
		[3,4,2,2,S(0.5)],
	]
	finalA = [zero(S),zero(S),zero(S),S(0.6)]
	
	B = [
		[1,2,2,3,S(0.3)],
		[2,3,1,2,S(0.4)],
		[3,3,1,2,S(0.6)]
	]
	finalB =S[zero(S),zero(S),S(0.7)]
	
	LA = coo2lil(A)	
	DA = coo2dense(S, A, length(LA), length(symbols))
	
	LB = coo2lil(B)
	DB = coo2dense(S, B, length(LB), length(symbols))
	
	DC = composition(DA,DB, length(symbols))
	C = dense2coo(DC)
	LC = coo2lil(C);
end

# ╔═╡ ce623bef-6f01-483c-b262-54d60e7751dc
begin
	fstA = VectorFST(
		LA,
		1,
		finalA
	)
	draw(fstA, symbols=symbols)
end

# ╔═╡ f0af788f-adb9-440c-b6b2-d0c5d8ad0943
begin
	fstB = VectorFST(
		LB,
		1,
		finalB
	)
	draw(fstB, symbols=symbols)
end

# ╔═╡ d0bbfe73-ba31-4cf5-9fb1-59ee2d636faf
begin
	fstC = VectorFST(
		LC,
		0,
		vec(external_dot(S, finalA, finalB))
	)
	draw(fstC, symbols=symbols)
end

# ╔═╡ Cell order:
# ╠═0a906f46-1ee0-4814-860a-56e127b71593
# ╠═8b30c1e2-f00d-4d8c-b91a-9d7401a4628f
# ╠═7e5488a5-4de3-4543-91f5-09e8db4f5c36
# ╠═314ba47a-2159-42ab-b580-f55759c75169
# ╠═5723bafe-c3e1-4b8f-8d22-7f9d6292f7f1
# ╠═699a88dc-ec13-4a4a-bdc9-f554a4be5e61
# ╠═807cde31-89ca-4b45-984d-f80d01ace895
# ╠═08e20858-4f43-4614-9d64-fcf7a32c7fa8
# ╠═c7b97138-1674-11ee-1f7b-e9882a147ade
# ╠═ce623bef-6f01-483c-b262-54d60e7751dc
# ╠═f0af788f-adb9-440c-b6b2-d0c5d8ad0943
# ╠═d0bbfe73-ba31-4cf5-9fb1-59ee2d636faf
