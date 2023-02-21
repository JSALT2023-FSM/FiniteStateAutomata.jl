# SPDX identifier: CECILL-2.1

"""
    abstract type AbstractFSA{K,L} end

Abstract base type for all finite state automata. `K` is the weight
semiring and `L` is the label type. Subtypes should implement the
following accessors.
"""
abstract type AbstractFSA{K,L} end

"""
    abstract type AbstractAcyclicFSA{K,L} <: AbstractFSA{K,L} end

Abstract base type for all cycle-free finite state automata.
"""
abstract type AbstractAcyclicFSA{K,L} <: AbstractFSA{K,L} end

"""
    semiring(A)

Return the semiring type of `A`.
"""
semiring(A::AbstractFSA{K}) where K = K

"""
    α(A)

Return the vector of initial states of `A`.
"""
α(::AbstractFSA)

"""
    T(A)

Return the transition matrix of `A`.
"""
T(::AbstractFSA)

"""
    ω(A)

Return the vector of final states of `A`.
"""
ω(::AbstractFSA)

"""
    ρ(A)

Return the weight of the emtpy string.
"""
ρ(::AbstractFSA)

"""
    λ(A)

Return the states' label of `A`.
"""
λ(::AbstractFSA)

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

"""
    I, V = initstates(A)

Return the indices of the initial states `I` and their associated value
`V`.
"""
initstates(A) = findnz(α(A))

"""
    I, J, V = finalstates(A)

Return the edges where `I` is the array of source states, `J` is the
array of destination states and `V` is the array of weights.
"""
edges(A::AbstractFSA) = findnz(T(A))

"""
    I, V = finalstates(A)

Return the indices of the final states `I` and their associated value
`V`.
"""
finalstates(A) = findnz(ω(A))

"""
    emptystring(A)

Return the weight of the empty string "ϵ" in `A`.
"""
emptystring(A) = ρ(A)


Base.parent(A::AbstractFSA) = A

"""
    convert(f::Function, A::AbstractFSA)

Convert the weight of `A`. The function `f` takes two arguments as
input: weight and the label of an edge and returnes a new semiring
value.
"""
Base.convert(f::Function, A::AbstractFSA)

#= Graph representation of a FSA using GraphViz =#

function Base.show(io::IO, ::MIME"image/svg+xml", A::AbstractFSA)
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

"""
    dot_write(io::IO, A::AbstractFSA)

Write a description of `A` on `io` in the [DOT](https://graphviz.org/doc/info/lang.html)
language.
"""
function dot_write(io::IO, A::AbstractFSA)
    println(io, "Digraph {")
    println(io, "rankdir=LR;")
    dot_write_nodes(io, T(A), λ(A))
    dot_write_initedges(io, α(A))
    dot_write_edges(io, T(A), ρ(A))
    dot_write_finaledges(io, ω(A))
    println(io, "}")
end

function dot_write_nodes(io::IO, T, λ)
    N = size(T, 1) # number of states
    println(io, "n0 [shape=\"point\"];")
    for i in 1:N
        #print(io, "n$(i) [label=\"$i|", λ[i], "\", shape=\"circle\"")
        print(io, "n$(i) [label=\"", escape_string("$(λ[i])"), "\", shape=\"circle\"")
        println(io, "];")
    end
    println(io, "n$(N+1) [shape=\"point\"];")
end

function dot_write_initedges(io::IO, α)
    for i in 1:size(α, 1)
        if ! iszero(α[i])
            dot_write_edge(io, 0, i, α[i])
        end
    end
end

function dot_write_edges(io, T, ρ)
    for i in 1:size(T, 1)
        for j in 1:size(T, 2)
            if ! iszero(T[i,j])
                dot_write_edge(io, i, j, T[i, j])
            end
        end
    end

    if ! iszero(ρ)
        dot_write_edge(io, 0, size(T, 2) + 1, ρ)
    end
end

function dot_write_finaledges(io::IO, ω)
    N = length(ω)
    for i in 1:size(ω, 1)
        if ! iszero(ω[i])
            dot_write_edge(io, i, N+1, ω[i])
        end
    end
end

function dot_write_edge(io, src, dst, weight)
    val = typeof(weight) <: AbstractFloat ? round(weight, digits=3) : weight
    print(io, "n$src -> n$dst [label=\"", weight, "\"];")
end

