# SPDX-License-Identifier: CECILL-2.1

struct DrawableFST{K,L} <: AbstractFST{K,L}
    fst::AbstractFST{K,L}
    isymbols::Dict
    osymbols::Dict
end

draw(fst::AbstractFST; isymbols = Dict(), osymbols = Dict()) =
    DrawableFST(fst, isymbols, osymbols)

function _show_svg(io::IO, fst::DrawableFST)
    dotpath, dotfile = mktemp()
    try
        fileio = IOContext(dotfile, :compact => true, :limit => true)
        show(fileio, MIME("text/vnd.graphviz"), fst)
        close(dotfile)
        svgpath, svgfile = mktemp()
        try
            run(`dot -Tsvg $(dotpath) -o $(svgpath)`)
            xml = read(svgfile, String)
            write(io, xml)
        finally
            close(svgfile)
            rm(svgpath)
        end
    finally
        rm(dotpath)
    end
end

Base.show(io::IO, ::MIME"image/svg+xml", fst::DrawableFST) = _show_svg(io, fst)

function Base.show(io::IO, ::MIME"text/vnd.graphviz", fst::DrawableFST)
    A = fst.fst
    isyms = fst.isymbols
    osyms = fst.osymbols

    println(io, "Digraph {")
    println(io, "rankdir=LR;")

    Q = nstates(A)
    P = narcs(A)

    for q in 1:nstates(A)
        style = "solid"
        shape = "circle"

        if ! iszero(α(A)[q])
            style = "bold"
        end

        if ! iszero(ω(A)[q])
            shape = "doublecircle"
        end

        if ! iszero(α(A)[q]) && ! isone(α(A)[q]) && ! iszero(ω(A)[q]) && ! isone(ω(A)[q])
            println(io, q, " [label=\"$q/", α(A)[q], "/", ω(A)[q], "\", style=\"$style\", shape=\"$shape\"];")
        elseif ! iszero(α(A)[q]) && ! isone(α(A)[q])
            println(io, q, " [label=\"$q/", A.α[q], "\", style=\"$style\", shape=\"$shape\"];")
        elseif ! iszero(ω(A)[q]) && ! isone(ω(A)[q])
            println(io, q, " [label=\"$q/", ω(A)[q], "\", style=\"$style\", shape=\"$shape\"];")
        else
            println(io, q, " [label=\"$q\", style=\"$style\", shape=\"$shape\"];")
        end
    end

    _get = (sym, l) -> get(sym, l, l)
    I, E1, V = findnz(S(A))
    E2, J, _ = findnz(A.D)
    p1 = sortperm(E1)
    p2 = sortperm(E2)
    for (s, d, w, a) in zip(I[p1], J[p2], V[p1], E1[p1])
        l = λ(A)[a]
        if l isa Pair
            label = join([_get(isyms, first(l)), _get(osyms, last(l))], ":")
        else
            label = _get(isyms, first(l))
        end

        if ! isone(w)
            println(io, s, "->", d, " [label=\"", label, "/", w, "\"];")
        else
            println(io, s, "->", d, " [label=\"", label, "\"];")
        end
    end

    println(io, "}")
end

