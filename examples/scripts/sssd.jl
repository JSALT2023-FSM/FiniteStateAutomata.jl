using ImageInTerminal
using Sixel
using FiniteStateAutomata

using NPZ



# Declare the semiring we want to use.
S = LogSemiring{Float32,1}

symbol_mapping = open(loadsymbols, "../assets/libri_examples/ctc_map_int.txt")
ctc_logits = npzread("../assets/libri_examples/2830-3980-0002/logits.npy")

Q, L = size(ctc_logits)
Q = Q + 1

W = S.(ctc_logits)

α = zeros(S, Q)
α[1] = one(1)

ω = zeros(S, Q)
ω[end] = one(S)

M = zeros(S, Q, Q, L, L)

for q in 1:Q-1
    for l in 1:L
        M[q, q+1, l, l] = W[q, l]
    end
end

ctc_fst = TensorFST(M, α, ω, (1, 2, 3, 4))

open("ctc.svg", "w") do f
    write(f, draw(ctc_fst, isymbols=symbol_mapping, osymbols=symbol_mapping) |> dot(:svg))
end


#This is for dense vectors
function dense_dot(x::Vector{S}, y::Vector{S}) where S
	semi_prod = x .⊗ y
	res = semi_prod[1]
	for i in 2:size(x)[1]
		res = res .⊕ semi_prod[i]
	end
	res
end


function permute_dim(A::FiniteStateAutomata.TensorFST)
	permuted_fst_matrix = deepcopy(A.M)	
	permuted_fst_matrix = permutedims(permuted_fst_matrix, [3, 4, 1, 2])

	return TensorFST(permuted_fst_matrix, A.α, A.ω, (3, 4, 1, 2))
	
end

function shortest_distance_fwd(A::FiniteStateAutomata.TensorFST)
	size_Q = numstates(A)
	T = reshape(reduce(hcat, sum(A.M, dims=[3;4])), size_Q, size_Q)
			
	u = A.α      # start vector
	SD = Matrix{S}[]
	#SD = vcat(SD, [u])
	#ones(S, numstates(ctc_fst), numstates(ctc_fst))
	#ρ = u .⊗ ctc_fst.ω ##this is ρ??
	#w = dense_dot(u, A.ω)
	

	while (length(findall(!iszero, u)) > 0)
		u = transpose(transpose(u) * T) 
		SD = vcat(SD, [u])				
	end
	sum(SD, dims=[1])[1]
end

function reweight(A::FiniteStateAutomata.TensorFST, is_fwd::Bool)
	Q = numstates(A)
	L = 
	if is_fwd
		SD = shortest_distance_fwd(deepcopy(A))
	else
		SD = shortest_distance_bwd(deepcopy(A))
	end

	diagonalized_sd = diagm(SD)
	diagonalized_sd_prev = diagm([inv(ele) for ele in vcat(SD[2:end], S[0.0])])
	

	i_l_len = size(A.M)[3]
	o_l_len = size(A.M)[4]
	reweighted_M = zeros(S, i_l_len, o_l_len, Q, Q)
	permuted_fst = permute_dim(A)

	#println(diagonalized_sd)
	#println(diagonalized_sd_prev)
	b = 0
	
	
	for n in 1:numstates(A)
		for m in 1:numstates(A)
			if is_fwd
				#reweighted_M[:, :, n, m] = diagonalized_sd * permuted_fst[:, :, n, m]
				#reweighted_M[:, :, n, m] = diagonalized_sd_prev * reweighted_M[:, :, n, m]
				b = diagonalized_sd * permuted_fst.M[:, :, n, m]
				b = diagonalized_sd_prev * b
				reweighted_M[:, :, n, m] = b
				
			end
		end
	end	
	reweighted_M
end

function shortest_distance_bwd(A::FiniteStateAutomata.TensorFST)
	size_Q = numstates(A)
	T′ = transpose(reshape(reduce(hcat, sum(A.M, dims=[3;4])), size_Q, size_Q))
			
	u = A.ω      # start vector
	SD = Matrix{S}[]
	#SD = vcat(SD, [u])
	#ones(S, numstates(ctc_fst), numstates(ctc_fst))
	#ρ = u .⊗ ctc_fst.ω ##this is ρ??
	#w = dense_dot(u, A.ω)
	

	while (length(findall(!iszero, u)) > 0)
		u = transpose(transpose(u) * T′) 
		SD = vcat(SD, [u])				
	end
	sum(SD, dims=[1])[1]
end

function transducer_weight(A::FiniteStateAutomata.TensorFST)
	#this uses the shortest distance but via the dot product gets to a number
	# that represents the weight up to that point
	sum_ = sum(A.M, dims=[3;4])
	size_Q = numstates(A)

	T = reshape(reduce(hcat, sum(A.M, dims=[3;4])), size_Q, size_Q)
			
	u = A.α      # start vector
	#ρ = u .⊗ ctc_fst.ω ##this is ρ??
	w = dense_dot(u, A.ω)

	while (length(findall(!iszero, u)) > 0)	
		u = transpose(transpose(u) * T)		
		w += dense_dot(u, A.ω)
	end
	w
end


function shortest_distance(A::FiniteStateAutomata.TensorFST)
	size_Q = numstates(A)
	T = reshape(reduce(hcat, sum(A.M, dims=[3;4])), size_Q, size_Q)
	u = A.α      # start vector
	SD = Matrix{S}[]

	while (length(findall(!iszero, u)) > 0)
		u = transpose(transpose(u) * T) 
		SD = vcat(SD, [u])				
	end
	SD
end




