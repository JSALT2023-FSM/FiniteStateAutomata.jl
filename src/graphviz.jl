# SPDX-License-Identifier: CECILL-2.1

struct DrawableFST{S,L}
    fst::AbstractFST{S,L}
    isymbols::Dict
    osymbols::Dict
    openfst_compat::Bool
end

draw(fst::AbstractFST; symbols = Dict(), isymbols = symbols, osymbols = symbols, openfst_compat = false) =
    DrawableFST(fst, isymbols, osymbols, openfst_compat)

_getlabel(l, isymbols, osymbols) = get(isymbols, l, l)
_getlabel(l::NTuple{1}, isymbols, osymbols) =  _getlabel(l[1], isymbols, osymbols)
_getlabel(l::NTuple{2}, isymbols, osymbols) = (
    _getlabel(l[1], isymbols, osymbols),
    _getlabel(l[2], isymbols, osymbols),
)
_getlabel(l::NTuple, isymbols, osymbols) = (
    _getlabel(l[1], isymbols, osymbols),
    _getlabel(l[2], isymbols, osymbols)...,
)
_getlabel(l::Pair, isymbols, osymbols) = join([
    get(isymbols, first(l), first(l)),
    get(osymbols, last(l), last(l)),
], ":")

function Base.show(io::IO, dfst::DrawableFST)
    fst = dfst.fst
    isyms = dfst.isymbols
    osyms = dfst.osymbols

    println(io, "Digraph {")
    println(io, "rankdir=LR;")

    offset = dfst.openfst_compat ? -1 : 0
    Q = nstates(fst)
    P = narcs(fst)

    for (q, iw, fw) in states(fst)
        q += offset
        style = "solid"
        shape = "circle"

        if ! iszero(iw)
            style = "bold"
        end

        if ! iszero(fw)
            shape = "doublecircle"
        end

        if ! iszero(iw) && ! isone(iw) && ! iszero(fw) && ! isone(fw)
            println(io, q, " [label=\"$q/", iw, "/", fw, "\", style=\"$style\", shape=\"$shape\"];")
        elseif ! iszero(iw) && ! isone(iw)
            println(io, q, " [label=\"$q/", iw, "\", style=\"$style\", shape=\"$shape\"];")
        elseif ! iszero(fw) && ! isone(fw)
            println(io, q, " [label=\"$q/", fw, "\", style=\"$style\", shape=\"$shape\"];")
        else
            println(io, q, " [label=\"$q\", style=\"$style\", shape=\"$shape\"];")
        end
    end

    for (s, d, l, w) in arcs(fst)
        s += offset
        d += offset

        label = escape_string(replace("$(_getlabel(l, isyms, osyms))", "\"" => ""))
        if ! isone(w)
            println(io, s, "->", d, " [label=\"", label, "/", w, "\"];")
        else
            println(io, s, "->", d, " [label=\"", label, "\"];")
        end
    end

    println(io, "}")
end


function _show_svg(io::IO, fst::DrawableFST)
    dotpath, dotfile = mktemp()
    try
        fileio = IOContext(dotfile, :compact => true, :limit => true)
        show(fileio, MIME("text/plain"), fst)
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

