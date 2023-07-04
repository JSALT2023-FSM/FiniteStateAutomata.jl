# SPDX-License-Identifier: CECILL-2.1

Base.convert(::Type{TensorFST{S, T}}, vfst::VectorFST{S}) where {S, T<:Array{S, 4}} = densefst(vfst)


# helper convertion functions 
# TODO: remove if unnceseary

# TODO: S should be read from A
function coo2arcs(A, S)
    t = Vector{Tuple{Int,Int,Int,S}}
    nstates = A.m
    states = Vector{t}()
    for i in 1:nstates
        push!(states,Vector{t}())
    end

    for (s,d,m) in zip(A.rows, A.cols, A.vals)
        for (i,j,v) in zip(m.rows, m.cols, m.vals)
            push!(states[s],(d,i,j,v))
        end
    end
    states
end

function vector2dict(A)
    # VectorFST to DictFST
    S = semiring(A)
    d = Dict{Int,Dict}()
    for (s,arcs) in enumerate(A.arcs)
        for arc in arcs
            if !haskey(d,s)
                d[s] = Dict{Int,Vector{Tuple{Int,Int,S}}}()
            end
            if !haskey(d[s], arc[1])
                d[s][arc[1]] = Vector{Tuple{Int,Int,S}}()
            end
            push!(d[s][arc[1]], arc[2:4])
        end
    end
    d
end

# TODO: dict_fst should be a proper struc where 
# we can get nstates and nsyms and the semiring S
function dict2coo(dict_fst, nstates, nsyms, S)
    M = SparseMatrixCOO{S,Int}
    state_rows = Vector{Int}()
    state_cols = Vector{Int}()
    state_vals = Vector{M}()
    for (s, sdic) in dict_fst
        for (d, ddic) in sdic
            push!(state_rows,s)
            push!(state_cols,d)
            label_rows = Vector{Int}()
            label_cols = Vector{Int}()
            label_vals = Vector{S}()
            for (i,j,v) in ddic
                push!(label_rows, i)
                push!(label_cols, j)
                push!(label_vals, v)
            end
            push!(state_vals, M(nsyms, nsyms, label_rows, label_cols, label_vals))
        end
    end
    SparseMatrixCOO{M, Int}(nstates, nstates, state_rows, state_cols, state_vals)
end
