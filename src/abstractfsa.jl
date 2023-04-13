# SPDX identifier: CECILL-2.1


function dot_write(io::IO, A::AbstractFSA)
    println(io, "Digraph {")
    println(io, "rankdir=LR;")
    dot_write_nodes(io, T(A), ω(A), ρ(A))
    dot_write_initedges(io, α(A), λ(A))
    dot_write_edges(io, T(A), ρ(A), λ(A))
    println(io, "}")
end

function dot_write_initedges(io::IO, α, λ)
    for i in 1:size(α, 1)
        if ! iszero(α[i])
            dot_write_edge(io, 0, i, λ[i], α[i])
        end
    end
end

function dot_write_edges(io, T, ρ, λ)
    for i in 1:size(T, 1)
        for j in 1:size(T, 2)
            if ! iszero(T[i,j])
                dot_write_edge(io, i, j, λ[j], T[i, j])
            end
        end
    end
end

function dot_write_edge(io, src, dst, label, weight)
    val = typeof(weight) <: AbstractFloat ? round(weight, digits=3) : weight
    print(io, "$src -> $dst [label=\"", label, "/", weight, "\"];")
end

