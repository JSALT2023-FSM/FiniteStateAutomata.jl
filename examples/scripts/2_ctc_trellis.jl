# SPDX-License-Identifier: CECILL-2.1

#=
Create a FST based on the output of a neural network trained with
the CTC-criterion for on utterance.
=#

using FiniteStateAutomata
using NPZ

# Load the symbol mapping and the logits.
symbol_mapping = open(loadsymbols, "../assets/libri_examples/ctc_map_int.txt")
ctc_logits = npzread("../assets/libri_examples/2830-3980-0002/logits.npy")

# Semiring to use. Possible choices are:
#   - LogSemiring{Float32|Float64,a} # a is the temperature parameter
#   - TropicalSemiring{Float32|Float64} # equivalent to LogSemiring{Float32|Float64,Inf}
S = LogSemiring{Float32,1}


# Number of states (# frames + 1) and the size of the alphabet (# tokens)[jk
Q, L = size(ctc_logits, 1)+1, size(ctc_logits, 2)

#=
NOTE: In the following, we use a very naive approach by allocating
dense tensor. You don't want to do that in practice as you will run
out of memory even for small FSTs.
=#

# Arcs of the FST stored in a 4-dimensional tensor with dimension :
# - 1: source state
# - 2: destination state
# - 3: input label
# - 4: output label
M = zeros(S, Q, Q, L, L)
for q in 1:Q-1
    for l in 1:L
        M[q, q+1, l, l] = S(ctc_logits[q, l])
    end
end

# Vector representing the starting state (only one initial state allowed).
α = zeros(S, Q)
α[1] = one(1)

# Vector of final weights.
ω = zeros(S, Q)
ω[end] = one(S)

# Create the TensorFST
ctc_fst = TensorFST(M, α, ω)

# Save the FST to a svg file.
path = "./ctc.svg"
println("Saving fst to $path")
open(path, "w") do f
    write(f, draw(ctc_fst, isymbols=symbol_mapping, osymbols=symbol_mapping) |> dot(:svg))
end

