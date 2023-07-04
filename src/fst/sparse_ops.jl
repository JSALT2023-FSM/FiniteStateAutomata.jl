# SPDX-License-Identifier: CECILL-2.1

#TODO remove nsymbols, missing getting number of symbols from vectorFST
#TODO when using A and B in tensor format conversion should not be necesary
function sparse_composition_kron(A, B, nsymbols)
	S = semiring(A)
	cooA = dict2coo(vector2dict(A), numstates(A),  nsymbols, S)
	cooB = dict2coo(vector2dict(B), numstates(B),  nsymbols, S)
	
	cooC = kron(cooA, cooB)
	
	arcsC = coo2arcs(cooC, S)
	initialA = zeros(S,numstates(A))
	initialB = zeros(S,numstates(B))
	initialA[A.initstate] = 1
	initialB[B.initstate] = 1
	initial = findfirst(x->x!=zero(S),kron(initialA, initialB))
	# TODO initial state
	VectorFST(arcsC, initial, kron(A.finalweights, B.finalweights) )
end