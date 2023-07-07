# SPDX-License-Identifier: CECILL-2.1

Base.convert(::Type{TensorFST{S, T}}, vfst::VectorFST{S}) where {S, T<:Array{S, 4}} = densefst(vfst)


# helper convertion functions 
# TODO: remove if unnceseary
# TODO: S should be read from A
function coo_sod2arcs(A, S)
    t = Vector{Tuple{Int,Int,Int,S}}
    nstates = A.m
    states = Vector{t}()
    for i in 1:nstates
        push!(states,Vector{t}())
    end

    for (s,d,m) in zip(A.rows, A.cols, A.vals)
        for (i,o,v) in zip(m.rows, m.cols, m.vals)
            push!(states[s],(d,i,o,v))
        end
    end
    states
end

function coo_lod2arcs(A, nstates, S)
    t = Vector{Tuple{Int,Int,Int,S}}  
    states = Vector{t}()
    for i in 1:nstates
        push!(states,Vector{t}())
    end

    for (i,o,m) in zip(A.rows, A.cols, A.vals)
        for (s,d,v) in zip(m.rows, m.cols, m.vals)
            push!(states[s],(d,i,o,v))
        end
    end
    states
end


# TODO: dict_fst should be a proper struc where 
# we can get nstates and nsyms and the semiring S
function dict2coo(dict_fst, outer_size, inner_size, S)
    M = SparseMatrixCOO{S,Int}
    outer_rows = Vector{Int}()
    outer_cols = Vector{Int}()
    outer_vals = Vector{M}()
    for (o1, o1dic) in dict_fst
        for (o2, o2dic) in o1dic
            push!(outer_rows,o1)
            push!(outer_cols,o2)
            inner_rows = Vector{Int}()
            inner_cols = Vector{Int}()
            inner_vals = Vector{S}()
            for (i1,i2,v) in o2dic
                push!(inner_rows, i1)
                push!(inner_cols, i2)
                push!(inner_vals, v)
            end
            push!(outer_vals, M(inner_size, inner_size, inner_rows, inner_cols, inner_vals))
        end
    end
    SparseMatrixCOO{M, Int}(outer_size, outer_size, outer_rows, outer_cols, outer_vals)
end

# VectorFST to DictFST_lod
function vector2dict_lod(A)
    S = semiring(A)
    d = Dict{Int,Dict}()
    for (s,arcs) in enumerate(A.arcs)
        for arc in arcs
			is = arc[2]
			os = arc[3]
            if !haskey(d, is)
                d[is] = Dict{Int,Vector{Tuple{Int,Int,S}}}()
            end
            if !haskey(d[is], os)
                d[is][os] = Vector{Tuple{Int,Int,S}}()
            end
            push!(d[is][os], (s, arc[1], arc[4]))
        end
    end
    d
end

# VectorFST to DictFST_sod
function vector2dict_sod(A)
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