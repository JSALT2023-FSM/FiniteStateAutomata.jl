# SPDX-License-Identifier: CECILL-2.1

dense_composition(A, B)	 = dense_composition_kron(A, B)

function dense_composition_matmul(A,B)
	S = semiring(A)
	MA, MB = M(A), M(B)

	n = numstates(A)
	m = numstates(B)

	ml = max(size(MA, 4), size(MB, 3))
	
	MA = cat(MA, 
		zeros(S, n, n, size(MA, 3), ml - size(MA, 4)); dims=4)
	MB = cat(MB, 
		zeros(S, m, m, ml - size(MB, 3), size(MB, 4) ); dims=3)
	MC = zeros(S,(n*m,n*m, size(MA,3), size(MB,4)))
	for i in 1:n
		for j in 1:n
			for k in 1:m
				for l in 1:m
					MC[(k-1)*n+i,(l-1)*n+j,:,:] .= MA[i,j,:,:]*MB[k,l,:,:]
				end
			end
		end
	end
	TensorFST(MC, kron(α(A), α(B)), kron(ω(A), ω(B)))
end

function dense_composition_kron(A, B)	
	S = semiring(A)
	MA, MB = M(A), M(B)

	pMA = permutedims(MA, (3, 4, 1, 2))
	pMB = permutedims(MB, (3, 4, 1, 2))

	ml = max(size(pMA, 2), size(pMB, 1))
	pMA = cat(pMA, 
		zeros(S, size(pMA, 1), ml - size(pMA, 2),numstates(A), numstates(A)); dims=2)
	pMB = cat(pMB, 
		zeros(S, ml - size(pMB, 1), size(pMB, 2), numstates(B), numstates(B)); dims=1)
	
	Q = numstates(A) * numstates(B)
	MC = zeros(S, size(pMA, 1), size(pMB, 2), Q, Q)
	for x in 1:size(MC, 1)
		for z in 1:size(MC, 2)
				MC[x,z,:,:] = sum(
					y -> kron(pMA[x,y,:,:], pMB[y,z,:,:]), 
					1:size(pMA, 2) ; 
					init=zeros(S, Q, Q)
				)
		end
	end

	TensorFST(permutedims(MC, (3, 4, 1, 2)), kron(α(A), α(B)), kron(ω(A), ω(B)))
end