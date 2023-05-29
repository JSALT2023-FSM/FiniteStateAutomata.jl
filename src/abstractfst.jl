# SPDX identifier: CECILL-2.1

const Label = Union{<:Integer,<:AbstractString}
const LabelMapping = Pair{<:Label,<:Label}

"""
    abstract type AbstractFST{K,L} end

Abstract base type for all WFST. `K` is the weight
semiring and `L` is the label type.
"""
abstract type AbstractFST{K,L} end

const Transducer = AbstractFST{<:Semiring,<:LabelMapping}
const Acceptor = AbstractFST{<:Semiring,<:Label}
const TransducerOrAcceptor = Union{<:Transducer,<:Acceptor}

"""
    α(A)

Return the vector of initial states of `A`.
"""
α(::TransducerOrAcceptor)

"""
    T(A)

Return the transition matrix of `A`.
"""
T(::TransducerOrAcceptor)

"""
    ω(A)

Return the vector of final states of `A`.
"""
ω(::TransducerOrAcceptor)

"""
    ρ(A)

Return the weight of the emtpy string.
"""
ρ(::TransducerOrAcceptor)

"""
    λ(A)

Return the states' label of `A`.
"""
λ(::TransducerOrAcceptor)

"""
    nstates(A)

Return the number of states in `A`.
"""
nstates(A) = length(λ(A))

"""
    nedges(A)

Return the number of edges in `A`.
"""
nedges(A) = nnz(T(A))

#= Graph representation of a FST using GraphViz =#

function Base.show(io::IO, ::MIME"image/svg+xml", A::AbstractFST)
    dotpath, dotfile = mktemp()
    try
        dot_write(IOContext(dotfile, :compact => true, :limit => true), A)
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

function Base.show(io::IO, ::MIME"image/svg+xml", tup::Tuple{<:AbstractFST, <:AbstractArray})
    A, highlights = tup
    dotpath, dotfile = mktemp()
    try
        dot_write(IOContext(dotfile, :compact => true, :limit => true), A; highlights)
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
    dot_write(io::IO, A::AbstractFST)

Write a description of `A` on `io` in the [DOT](https://graphviz.org/doc/info/lang.html)
language.
"""
function dot_write(io::IO, A::AbstractFST; highlights = [], hcolor = "blue")
    println(io, "Digraph {")
    println(io, "rankdir=LR;")
    dot_write_nodes(io, T(A), ω(A), ρ(A))
    dot_write_initedges(io, α(A), λ(A))
    dot_write_transmat(io, T(A), λ(A))
    #dot_write_edges(io, T(A), ρ(A), λ(A))

    if ! isempty(highlights)
        println(io, join(highlights, ","), " [style=\"bold\", color=\"$hcolor\"];")
    end
    println(io, "}")
end

function dot_write_transmat(io::IO, A::TransitionMatrix, λ::AbstractVector)
    Q = size(A.S, 1)

    I, J, V = findnz(A.S)
    for (i, j, v) in zip(I, J, V)
        _dot_write_edge(io, i, j, λ[j], v)
    end

    I, J, V = findnz(A.U)
    for (i, j, v) in zip(I, J, V)
        _dot_write_edge(io, i, Q+j, "ϵ", v)
    end

    I, J, V = findnz(A.E)
    for (i, j, v) in zip(I, J, V)
        _dot_write_edge(io, Q+i, Q+j, "ϵ", v)
    end

    I, J, V = findnz(A.V)
    for (i, j, v) in zip(I, J, V)
        _dot_write_edge(io, Q+i, j, λ[j], v)
    end
end

function dot_write_edge(io, src, dst, label, weight)
    if isone(weight)
        println(io, "$src -> $dst [label=\"", label, "\"];")
    else
        style = iszero(weight) ? "style=invis" : ""
        val = typeof(weight) <: AbstractFloat ? round(weight, digits=3) : weight
        println(io, "$src -> $dst [label=\"", label, "/", weight, "\"", style, "];")
    end
end

function dot_write_edge(io, src, dst, label::Pair, weight)
    if isone(weight)
        println(io, "$src -> $dst [label=\"", label[1], ":", label[2], "\"];")
    else
        style = iszero(weight) ? "style=invis" : ""
        val = typeof(weight) <: AbstractFloat ? round(weight, digits=3) : weight
        println(io, "$src -> $dst [label=\"", label[1], ":", label[2], "/", weight, "\" ", style, "];")
    end
end

function dot_write_nodes(io::IO, T, ω, ρ)
    Q, P = size(T.S, 1), size(T.E, 1)
    println(io, join(0:Q+P, ","), " [shape=\"circle\"];")

    if iszero(ρ)
        println(io, "0 [style=\"bold\"];")
    else
        println(io, "0 [label=\"0/", ρ, "\", shape=\"doublecircle\", style=\"bold\"];")
    end

    for i in 1:size(T, 1)
        if ! iszero(ω[i])
            if ! isone(ω[i])
                println(io, "$(i) [label=\"$i/", ω[i], "\", shape=\"doublecircle\"];")
            else
                println(io, "$(i) [label=\"$i\", shape=\"doublecircle\"];")
            end
        end
    end
end

function dot_write_initedges(io::IO, α, λ)
    for i in 1:size(α, 1)
        if ! iszero(α[i])
            dot_write_edge(io, 0, i, λ[i], α[i])
        end
    end
end

function dot_write_edges(io, T::MatrixPowerSum, ρ, λ)
    I, J, V = findnz(T)
    for (i, j, w) in zip(I, J, V)
        dot_write_edge(io, i, j, λ[j], w)
    end
end

