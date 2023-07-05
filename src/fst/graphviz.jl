# SPDX-License-Identifier: CECILL-2.1

#=====================================================================#
# Write the FST in dot graphviz format.
#=====================================================================#

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

#=====================================================================#
# Backends for the dot graphviz representation:
#   - svg
#=====================================================================#

struct DrawableFST
    fst
    isyms
    osyms
    openfst_compat
end

function draw(fst::AbstractFST; symbols = Dict(), isymbols = symbols, osymbols = symbols, openfst_compat = false)
    DrawableFST(fst, isymbols, osymbols, openfst_compat)
end

function _show_svg(io, fst, isyms, osyms, openfst_compat)
    dotpath, dotfile = mktemp()
    fileio = IOContext(dotfile, :compact => true, :limit => true)
    _write_dot(fileio, fst, isyms, osyms, openfst_compat)
    close(dotfile)
    svgpath, svgfile = mktemp()
    run(`dot -Tsvg $(dotpath) -o $(svgpath)`)
    xml = read(svgfile, String)
    write(io, xml)
    rm(svgpath)
end

Base.show(io::IO, ::MIME"image/svg+xml", fst::DrawableFST) =
        _show_svg(io, fst.fst, fst.isyms, fst.osyms, fst.openfst_compat)

Base.show(io::IO, ::MIME"text/vnd.graphviz", fst::DrawableFST) =
    _write_dot(IOContext(io, :compact => true), fst.fst, fst.isyms, fst.osyms, fst.openfst_compat)

