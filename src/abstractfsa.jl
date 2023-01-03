
"""
    abstract type AbstractFSA{K,L} end

Abstract base type for all finite state automata. `K` is the weight
semiring and `L` is the label type.
"""
abstract type AbstractFSA{K,L} end

#= Graph representation of a FSA using GraphViz =#

function Base.show(io::IO, ::MIME"image/svg+xml", fsa::FSA)
    dotfsa = dot(fsa)
    dotpath, dotfile = mktemp()
    try
        write(dotfile, dotfsa)
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

"""
    dot(fsa)

Return a string describing `fsa` in the [DOT](https://graphviz.org/doc/info/lang.html)
language.
"""
function dot(fsa::AbstractFSA)
    out = "Digraph {\n"
    out *= "rankdir=LR;\n"
    out = dot_write_nodes(out, fsa.T, fsa.λ)
    out = dot_write_initedges(out, fsa.α)
    out = dot_write_edges(out, fsa.T, fsa.λ)
    out = dot_write_finaledges(out, fsa.ω)
    out *= "}\n"
end

function dot_write_nodes(out, T, λ)
    N = size(T, 1) # number of states
    out *= "n0 [shape=\"point\"];\n"
    for i in 1:N
        out *= "n$(i) [label=$(λ[i]), shape=\"circle\""
        out *= "];\n"
    end
    out *= "n$(N+1) [shape=\"point\"];\n"
end

function dot_write_initedges(out, α)
    for i in 1:size(α, 1)
        if ! iszero(α[i])
            out = dot_write_edge(out, 0, i, α[i])
        end
    end
    out
end

function dot_write_edges(out, T, λ)
    for i in 1:size(T, 1)
        for j in 1:size(T, 2)
            if ! iszero(T[i,j])
                out = dot_write_edge(out, i, j, T[i, j])
            end
        end
    end
    out
end

function dot_write_finaledges(out, ω)
    N = length(ω)
    for i in 1:size(ω, 1)
        if ! iszero(ω[i])
            out = dot_write_edge(out, i, N+1, ω[i])
        end
    end
    out
end

function dot_write_edge(out, src, dst, weight)
    out *= "n$src -> n$dst [label=\""
    #out *= isone(weight) ? "" : "/$(round(weight, digits=3))"
    out *= "$(round(weight, digits=3))"
    out *= "\"];\n"
end

