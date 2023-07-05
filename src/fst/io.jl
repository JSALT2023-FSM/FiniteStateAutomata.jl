# SPDX-License-Identifier: CECILL-2.1

"""
    draw(fst[; isymbols=Dict(), osymbols=Dict(), openfst_compat=false])

Return a string containg the dot formatted representation of `fst`.

For a quick display in Pluto notebook:

    draw(fst[; kwargs...]) |> dot([backend=svg|png]) |> HTML
"""
function draw(fst::AbstractFST;
              isymbols = Dict(),
              osymbols = Dict(),
              openfst_compat = false)

    buffer = IOBuffer()
    ctx = IOContext(buffer, :compact => true)
    _write_dot(ctx, fst, isymbols, osymbols, openfst_compat)
    String(take!(buffer))
end

_getlabel(il, ol, isymbols, osymbols) = join([
    get(isymbols, il, il),
    get(osymbols, ol, ol),
], ":")

function _write_dot(io::IO, fst, isyms, osyms, openfst_compat)
    println(io, "Digraph {")
    println(io, "rankdir=LR;")

    offset = openfst_compat ? -1 : 0

    for q in states(fst)
        fw = finalweight(fst, q)
        q += offset
        print(io, q + offset, " [label=\"")
        (iszero(fw) || isone(fw)) ? print(io, q) : print(io, q, "/", fw)
        print(io, "\", style=\"", isinit(fst, q) ? "bold" : "solid", "\", ")
        println(io, "shape=\"", iszero(fw) ? "circle" : "doublecircle", "\"];")
    end

    for s in states(fst)
        for (d, il, ol, w) in arcs(fst, s)
            print(io, s + offset, " -> ", d + offset, " [label=\"")
            print(io, get(isyms, il, il), ":", get(osyms, ol, ol))
            ! isone(w) ? println(io, "/", w, "\"];") : println(io, "\"];")
        end
    end

    println(io, "}")
end

struct Dot2SVG end
function (::Dot2SVG)(str)
    dir = mktempdir()
    dotpath = joinpath(dir, "fst.dot")
    svgpath = joinpath(dir, "fst.svg")
    open(f -> write(f, str), dotpath, "w")
    run(`dot -Tsvg $(dotpath) -o $(svgpath)`)
    open(f -> read(f, String), svgpath)
end

struct Dot2PNG end
function (::Dot2PNG)(str)
    dir = mktempdir()
    dotpath = joinpath(dir, "fst.dot")
    pngpath = joinpath(dir, "fst.png")
    open(f -> write(f, str), dotpath, "w")
    run(`dot -Tpng $(dotpath) -o $(pngpath)`)
    load(pngpath)
end

dot(backend) = dot(Val(backend))
dot(::Val{:svg}) = Dot2SVG()
dot(::Val{:png}) = Dot2PNG()


"""
    loadsymbols(str)

Load a symbol from a string `str` formatted as
```
ϵ  0
sym1 1
sym2 2
sym3 3
...
```

!!! note
    The 0 identifier is reserved for the special symbol ``\\epsilon``.

"""
function loadsymbols(lines)
	symtable = Dict()
	for line in lines
		tokens = split(line)
		length(tokens) == 0 && continue # skip emtpy lines
		label, id = tokens
		symtable[parse(Label, id)] = label
	end
	symtable
end

loadsymbols(str::AbstractString) = loadsymbols(split(str, "\n"))
loadsymbols(io::IOStream)= loadsymbols(eachline(io))

"""
    compile(wfst[; semiring = LogSemiring{Float32,1}, acceptor = false, openfst_compat = false])

Create a `SparseFST` object from a text-formatted FST definition file.
"""
function compile(lines;
                 semiring=LogSemiring{Float32,1},
                acceptor=false,
                openfst_compat=false)
    offset = openfst_compat ? 1 : 0

    S = semiring
    arcs = Dict{State,Vector{Arc{S}}}()
    fstates = State[]
    fweights = S[]
    maxstate = 0
    initstate = 0

    init = false
    for (i, line) in enumerate(lines)
        tokens = split(line)
        isempty(tokens) && continue

        if ! init
            initstate = parse(State, tokens[1]) + offset
            init = true
        end

        if 1 ≤ length(tokens) ≤ 2
            state = parse(State, tokens[1]) + offset
            weight = length(tokens) == 1 ? one(S) : S(parse(Float64, tokens[2]))
            push!(fstates, state)
            push!(fweights, weight)
            maxstate = max(maxstate, state)
        else
            src = parse(State, tokens[1]) + offset
            dest = parse(State, tokens[2]) + offset
            maxstate = max(maxstate, max(src,dest))
            isym = parse(Label, tokens[3])
            if acceptor
                osym = isym
                weight = length(tokens) == 3 ? one(S) : S(parse(Float64, tokens[4]))
            else
                osym = parse(State, tokens[4])
                weight = length(tokens) == 4 ? one(S) : S(parse(Float64, tokens[5]))
            end

            arcs[src] = push!(
                get(arcs, src, Arc{S}[]),
                (dest, isym, osym, weight)
            )
        end
    end


    arclist = [get(arcs, src, Arc{S}[]) for src in 1:maxstate]
    finalweights = zeros(S, maxstate)
    finalweights[fstates] .= fweights
    VectorFST(arclist, initstate, finalweights)
end
compile(io::IOStream; kwargs...) = compile(eachline(io); kwargs...)
compile(txt::AbstractString; kwargs...) = compile(split(txt, "\n"); kwargs...)

function Base.print(io::IO, fst::AbstractFST; openfst_compat = false, acceptor = false)
    offset = openfst_compat ? -1 : 0
    for s in states(fst)
        for (d, il, ol, w) in arcs(fst, s)
            s += offset
            d += offset
            if acceptor
                print(io, s, " ", d, " ", il)
                isone(w) ? println(io) : println(io, " ", w)
            else
                print(io, s, " ", d, " ", il, " ", ol)
                isone(w) ? println(io) : println(io, " ", w)
            end
        end
    end

    nfs = 0
    for s in states(fst)
        if ! iszero(finalweight(fst, s))
            nfs += 1
        end
    end

    count = 1
    for s in states(fst)
        fw = finalweight(fst, s)
        if ! iszero(fw)
            isone(fw) ? print(io, s) : print(io, s, " ", fw)
            if count < nfs
                println(io)
            end
            count += 1
        end
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