# SPDX-License-Identifier: CECILL-2.1

function _fsadet(M, edges, states)
	state2idx = Dict(q => i for (i, q) in enumerate(sort(collect(states))))
    D, Q = length(states), length(λ(M))
    K = eltype(α(M))

	I, V = [], K[]
	for i in edges[()]
		eᵢ = sparsevec(collect(i), one(K), Q)
		push!(I, state2idx[i])
        push!(V, dot(α(M), eᵢ))
	end
	a = sparsevec(I, V, D)

	I, J, V = [], [], K[]
	for i in keys(edges)
		i == () && continue
		for j in edges[i]
			push!(I, state2idx[i])
			push!(J, state2idx[j])
			eᵢ = sparsevec(collect(i), one(K), Q)
			eⱼ = sparsevec(collect(j), one(K), Q)
            push!(V, eᵢ' * T(M) * eⱼ)
		end
	end
	B = sparse(I, J, V, D, D)

	I, V = [], K[]
	for i in states
		eᵢ = sparsevec(collect(i), one(K), Q)
        v = dot(ω(M), eᵢ)
        if v ≠ zero(K)
            push!(I, state2idx[i])
            push!(V, v)
        end
	end
	o = sparsevec(I, V, D)

    l = [λ(M)[q[1]] for q in sort(collect(states))]

    FST(a, B, o, ρ(M), l)
end

function _spmap(K, mapping)
	unique_ids = Set(mapping)
	id2idx = Dict(id => i for (i, id) in enumerate(collect(unique_ids)))
	idx_mapping = [id2idx[id] for id in mapping]

	Q, D = length(mapping), length(unique_ids)
	sparse(1:Q, idx_mapping, one(K), Q, D)
end

"""
    determinize(A::AbstractFST)

Return an "equivalent" deterministic version of `A`. Note that the
weights of the merged paths are simply added and locally renormalized.
Therefore, the resulting FST may have different weighting than the
original FST.
"""
function determinize(A::AbstractFST{K}) where K
    C = _spmap(K, λ(A))

	# Edges of the determinized FST.
	edges = Dict()

	# Visited states of the determinized FST.
	visited = Set()

	# Initialize the stack, `()` represents the initial state
	# of the determinized FST.
    stack = Any[((), α(A))]
	while ! isempty(stack)
		# Get the first element from the stack.
		# `a` is the ancestor state (of the new FST)
		# and `v` is the children of `a` in `A`.
		a, v = pop!(stack)

		# Merge states with the same label.
		V = spdiagm(v) * C

		# Iterate over the column of V that have non-zero elements.
		for i in 1:size(V, 2)
			c = V[:, i]
			nnz(c) == 0 && continue

			# A new state is defined based on the indices
			# of the non-zero elements.
			q = Tuple(findnz(c)[1])

			# Record the edge between the ancestor and the new
			# state.
			edges[a] = push!(get(edges, a, []), q)

			# If the new state `q` has not been visited yet
			# add its children (represented as a vector)
			# on the stack.
			if q ∉ visited
                push!(stack, (q, T(A)' * c))
				push!(visited, q)
			end
		end
	end
    _fsadet(A, edges, visited)
end

"""
    minimize(A::AbstractFST)

Return a minimal equivalent FST.
"""
minimize(A::AbstractFST) = (reverse ∘ determinize ∘ reverse ∘ determinize)(A)

