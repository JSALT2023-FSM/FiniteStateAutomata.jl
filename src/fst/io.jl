# SPDX-License-Identifier: CECILL-2.1

"""
    symboltable(textdata)

Load a symbol table.
"""
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

"""
    compile(wfst[; semiring = LogSemiring{Float32,1}, acceptor = false, openfst_compat = false])

Create a `SparseFST` object from a text-formatted FST definition file.
"""
function compile(wfst::AbstractString; semiring = LogSemiring{Float32,1},
                acceptor = false, openfst_compat = false)
    offset = openfst_compat ? 1 : 0

    K = semiring
    arcs = acceptor ? Dict{Int,Any}() : Dict{Pair{Int, Int},Any}()
    fstates = Int[]
    fweights = K[]
    Qset = Set{Int}()

    lines = split(wfst, "\n")
    for line in lines
        tokens = split(line)
        isempty(tokens) && continue

        if 1 ≤ length(tokens) ≤ 2
            state = parse(Int, tokens[1]) + offset
            weight = length(tokens) == 1 ? one(K) : K(parse(Float64, tokens[2]))
            push!(fstates, state)
            push!(fweights, weight)
            push!(Qset, state)
        else
            src = parse(Int, tokens[1]) + offset
            dest = parse(Int, tokens[2]) + offset
            push!(Qset, src)
            push!(Qset, dest)

            isym = parse(Int, tokens[3])
            if ! acceptor
                osym = parse(Int, tokens[4])
                weight = length(tokens) == 4 ? one(K) : K(parse(Float64, tokens[5]))
                arcs[isym => osym] = push!(
                    get(arcs, isym => osym, []),
                    (src, dest, weight)
                )
            else
                weight = length(tokens) == 3 ? one(K) : K(parse(Float64, tokens[4]))
                arcs[isym] = push!(
                    get(arcs, isym => osym, []),
                    (src, dest, weight)
                )
            end
        end
    end

    initstate = parse(Int, split(first(lines))[1]) + offset
    Q = length(Qset)

    function spm(arclist, Q)
        I, J, V = [], [], K[]
        for (s, d, v) in arclist
            push!(I, s)
            push!(J, d)
            push!(V, v)
        end
        sparse(I, J, V, Q, Q)
    end

    λ = collect(keys(arcs))
    M = SparseMatrices([spm(arcs[l], Q) for l in λ]...)

    SparseFST(
        M,
        sparsevec([initstate], one(K), Q),
        sparsevec(fstates, fweights, Q),
        λ
    )
end

function Base.print(io::IO, fst::AbstractFST; openfst_compat = false)
    if nnz(α(fst)) != 1 && ! isone(nonzeros(α(fst))) > 1
        throw(ArgumentError("Can only print FST with a unique starting state with initial weight 1̄"))
    end

    offset = openfst_compat ? -1 : 0

    for (s, d, l, w) in arcs(fst)
        s += offset
        d += offset
        if fst isa Acceptor
            println(io, s, " ", d, " ", l, " ", w)
        else
            println(io, s, " ", d, " ", first(l), " ", last(l), " ", w)
        end
    end

    I, V = findnz(ω(fst))
    for (q, w) in zip(I, V)
        q += offset
        println(io, q, " ", w)
    end
end

Base.print(fst::AbstractFST; openfst_compat::Bool = false) =
    print(IOContext(stdout, :compact => true), fst; openfst_compat)

struct FSTSummary
    fst::AbstractFST
end

Base.summary(fst::AbstractFST) = FSTSummary(fst)

function Base.show(io::IO, ::MIME"text/html", s::FSTSummary)
    fst = s.fst
    println(io, "<table>")
    println(io, "<tr><td>type</td><td>", fst isa Acceptor ? "acceptor" : "transducer", "</td></tr>")
    println(io, "<tr><td>storage</td><td>", typeof(fst).name.name, "</td></tr>")
    println(io, "<tr><td>label type</td><td>", eltype(λ(fst)), "</td></tr>")
    println(io, "<tr><td>semiring</td><td>", semiring(fst), "</td></tr>")
    println(io, "<tr><td># arcs</td><td>", narcs(fst), "</td></tr>")
    println(io, "<tr><td># epsilon arcs</td><td>",
            length(findall(l -> isepsilon(l), λ(fst))),
            "</td></tr>")
    println(io, "<tr><td># input epsilon arcs</td><td>",
            length(findall(l -> isinputepsilon(l), λ(fst))),
            "</td></tr>")
    println(io, "<tr><td># output epsilon arcs</td><td>",
            length(findall(l -> isoutputepsilon(l), λ(fst))),
            "</td></tr>")
    println(io, "<tr><td># states</td><td>", nstates(fst), "</td></tr>")
    println(io, "<tr><td># accessible states</td><td>", sum(accessible(fst)), "</td></tr>")
    println(io, "<tr><td># coaccessible states</td><td>", sum(coaccessible(fst)), "</td></tr>")
    print(io, "</table>")
end

