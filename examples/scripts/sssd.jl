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

begin
	sum_ = sum(ctc_fst.M, dims=[3;4])
	size_Q = numstates(ctc_fst)

	T = reshape(reduce(hcat, sum(ctc_fst.M, dims=[3;4])), size_Q, size_Q)
			
	n = 0	
	u = ctc_fst.α      # start vector
	#ρ = u .⊗ ctc_fst.ω ##this is ρ??
	w = dense_dot(u, ctc_fst.ω)

	while (n < size_Q)
		n += 1		
		u = transpose(transpose(u) * T)		
		w += dense_dot(u, ctc_fst.ω)
		
	end

	w

end




