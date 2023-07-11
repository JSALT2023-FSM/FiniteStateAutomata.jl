using FiniteStateAutomata
using BenchmarkTools
using Profile
using ProfileSVG
using PProf
using Base.Threads

@show nthreads()

K = LogSemiring{Float16, ℯ}

begin
	label_mapping = open("examples/assets/libri_examples/ctc_map.txt") do f
		readlines(f)
	end
	label_mapping = Dict(k => v for (k,v) in enumerate(label_mapping))
	label_mapping_reverse = Dict(value => key for (key, value) in label_mapping)
	nsym = length(label_mapping)
end;

A = open(x-> compile(x,semiring=K, openfst_compat=false), "examples/assets/logits.fst.txt") 
B = open(x-> compile(x,semiring=K, openfst_compat=true), "examples/assets/charwotw.3.rmeps.fst.txt") 
Q = numstates(A) * numstates(B)

cooAsod = dict2coo_csc(vector2dict_sod(A), numstates(A), nsym, K)
cooBsod = dict2coo_csc(vector2dict_sod(B), numstates(B), nsym, K)

Profile.clear()
@pprof kron_coo_mt(cooAsod, cooBsod)

# display(@benchmark kron_coo_mt(cooAsod, cooBsod))

cooAlod = dict2coo(vector2dict_lod(A), nsym, numstates(A), K)
cooBlod = dict2coo(vector2dict_lod(B), nsym, numstates(B), K)

display(@benchmark sparse_coo_composition_lod_mt(cooAlod, cooBlod, K, nsym, Q))
