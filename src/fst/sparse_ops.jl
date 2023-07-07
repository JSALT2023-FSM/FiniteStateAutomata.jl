# SPDX-License-Identifier: CECILL-2.1

#TODO remove nsymbols, missing getting number of symbols from vectorFST
#TODO when using A and B in tensor format conversion should not be necesary

# sod goes for state outermost ordering
function sparse_composition_sod(A, B, nsymbols)
	S = semiring(A)
	cooA = dict2coo(vector2dict_sod(A), numstates(A),  nsymbols, S)
	cooB = dict2coo(vector2dict_sod(B), numstates(B),  nsymbols, S)
	
	cooC = kron(cooA, cooB)	
	arcsC = coo_sod2arcs(cooC, S)

	initialA = zeros(S,numstates(A))
	initialB = zeros(S,numstates(B))
	initialA[A.initstate] = 1
	initialB[B.initstate] = 1
	initial = findfirst(x->x!=zero(S),kron(initialA, initialB))
	# TODO initial state
	VectorFST(arcsC, initial, kron(A.finalweights, B.finalweights) )
end

# lod goes for labels in outermost dimension
function sparse_composition_lod(A, B, nsymbols)	
	S = semiring(A)
    M = SparseMatrixCOO{S,Int}	
	
	cooA = dict2coo(vector2dict_lod(A), nsymbols, numstates(A), S)
	cooB = dict2coo(vector2dict_lod(B), nsymbols, numstates(B), S)
	
	Q = numstates(A) * numstates(B)
	
    label_rows = Vector{Int}()
    label_cols = Vector{Int}()
    label_vals = Vector{M}()

	for x in 1:nsymbols
		for z in 1:nsymbols	
			# state_c = sparse([1],[1],[zero(S)],Q,Q)
			kronresult_list = []
			for y in 1:nsymbols
				# @show y, cooA[x,y], cooB[y,z], hasitem(cooA,x,y), hasitem(cooB,y,z)
				if hasitem(cooA,x,y) && hasitem(cooB,y,z)					
					kronresult = kron(cooA[x,y], cooB[y,z])
					# TODO quick version coo to csc to do addition, should handle proper coo+coo 
					push!(kronresult_list, sparse(kronresult.rows, kronresult.cols, kronresult.vals, Q,Q) )
				end
			end
			if length(kronresult_list)>0
				push!(label_rows,x)
				push!(label_cols,z)	
				push!(label_vals,tocoo(sum(kronresult_list)))
			end
		end		
	end

	cooC = SparseMatrixCOO{M, Int}(nsymbols, nsymbols, label_rows, label_cols, label_vals)
	arcsC = coo_lod2arcs(cooC, Q, S)

	initialA = zeros(S,numstates(A))
	initialB = zeros(S,numstates(B))
	initialA[A.initstate] = S(1)
	initialB[B.initstate] = S(1)
	initial = findfirst(x->x!=zero(S),kron(initialA, initialB))
	final = kron(A.finalweights, B.finalweights) 

	VectorFST(arcsC, initial, final )
	
end