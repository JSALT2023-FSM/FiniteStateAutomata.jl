# SPDX-License-Identifier: CECILL-2.1

function symboltable(txt::AbstractString)
	symtable = Dict()
	for line in split(txt, "\n")
		tokens = split(line)
		length(tokens) == 0 && continue # skip emtpy lines
		label, id = tokens
		symtable[parse(Int, id)] = label
	end
	symtable
end

function compile(wfst::AbstractString;
                 semiring = LogSemiring{Float32},
                 acceptor = false,
                 isymbols = Dict(),
                 osymbols = Dict())
    K = semiring
    arcslabel = acceptor ? Int[] : Pair{Int,Int}[]
    arcsweight = K[]
    arcssrc = Int[]
    arcsdest = Int[]
    fstates = Int[]
    fweights = K[]

    lines = split(wfst, "\n")
    for line in lines
        tokens = split(line)
        isempty(tokens) && continue

        if 1 ≤ length(tokens) ≤ 2
            state = parse(Int, tokens[1])
            weight = length(tokens) == 1 ? one(K) : K(parse(Float64, tokens[2]))
            push!(fstates, state)
            push!(fweights, weight)
        else
            src = parse(Int, tokens[1])
            dest = parse(Int, tokens[2])
            push!(arcssrc, src)
            push!(arcsdest, dest)

            isym = parse(Int, tokens[3])
            if ! acceptor
                osym = parse(Int, tokens[4])
                weight = length(tokens) == 4 ? one(K) : K(parse(Float64, tokens[5]))
                push!(arcslabel, isym => osym)
            else
                weight = length(tokens) == 3 ? one(K) : K(parse(Float64, tokens[4]))
                push!(arcslabel, isym)
            end
            push!(arcsweight, weight)
        end
    end

    Qset = Set{Int}(arcssrc) ∪ Set{Int}(arcsdest) ∪ Set{Int}(fstates)

    initstate = parse(Int, split(first(lines))[1])
    Q = length(Qset)
    A = length(arcslabel)

    FST(
        sparsevec([initstate], one(K), Q),
        sparse(arcssrc, collect(1:A), arcsweight, Q, A),
        sparse(collect(1:A), arcsdest, one(K), A, Q),
        sparsevec(fstates, fweights, Q),
        arcslabel,
        isymbols,
        osymbols
    )
end

function Base.print(io::IO, fst::AbstractFST)
    if nnz(α(fst)) != 1 && ! isone(nonzeros(α(fst))) > 1
        throw(ArgumentError("Can only print FST with a unique starting state with initial weight 1̄"))
    end

    I, A1, V = findnz(S(fst))
    A2, J, _ = findnz(D(fst))
    p1 = sortperm(A1)
    p2 = sortperm(A2)
    for (s, d, w, a) in zip(I[p1], J[p2], V[p1], A1[p1])
        if fst isa Acceptor
            println(io, s, " ", d, " ", λ(fst)[a], " ", w)
        else
            println(io, s, " ", d, " ", first(λ(fst)[a]), " ", last(λ(fst)[a]), " ", w)
        end
    end

    I, V = findnz(ω(fst))
    for (q, w) in zip(I, V)
        println(io, q, " ", w)
    end
end

