# # SPDX-License-Identifier: CECILL-2.1

# struct DrawableFST{S}
#     fst::AbstractFST{S}
#     isymbols::Dict
#     osymbols::Dict
#     openfst_compat::Bool
# end

# draw(fst::AbstractFST; symbols = Dict(), isymbols = symbols, osymbols = symbols, openfst_compat = false) =
#     DrawableFST(fst, isymbols, osymbols, openfst_compat)

# _getlabel(il, ol, isymbols, osymbols) = join([
#     get(isymbols, il, il),
#     get(osymbols, ol, ol),
# ], ":")

# function Base.show(io::IO, dfst::DrawableFST)
#     fst = dfst.fst
#     isyms = dfst.isymbols
#     osyms = dfst.osymbols

#     println(io, "Digraph {")
#     println(io, "rankdir=LR;")

#     offset = dfst.openfst_compat ? -1 : 0

#     for q in states(fst)
#         fw = finalweight(fst, q)
#         q += offset
#         print(io, q + offset, " [label=\"")
#         (iszero(fw) || isone(fw)) ? print(io, q) : print(io, q, "/", fw)
#         print(io, "\", style=\"", isinit(fst, q) ? "bold" : "solid", "\", ")
#         println(io, "shape=\"", iszero(fw) ? "circle" : "doublecircle", "\"];")
#     end

#     for s in states(fst)
#         for (d, il, ol, w) in arcs(fst, s)
#             print(io, s + offset, " -> ", d + offset, " [label=\"")
#             print(io, isyms[il], ":", osyms[ol])
#             ! isone(w) ? print(io, "/", w, "\"];") : print(io, "\"];")
#         end
#     end

#     println(io, "}")
# end

# function _show_svg(io::IO, fst::DrawableFST)
#     dotpath, dotfile = mktemp()
#     try
#         fileio = IOContext(dotfile, :compact => true, :limit => true)
#         show(fileio, MIME("text/plain"), fst)
#         close(dotfile)
#         svgpath, svgfile = mktemp()
#         try
#             run(`dot -Tsvg $(dotpath) -o $(svgpath)`)
#             xml = read(svgfile, String)
#             write(io, xml)
#         finally
#             close(svgfile)
#             rm(svgpath)
#         end
#     finally
#         rm(dotpath)
#     end
# end

# Base.show(io::IO, ::MIME"image/svg+xml", fst::DrawableFST) = _show_svg(io, fst)

