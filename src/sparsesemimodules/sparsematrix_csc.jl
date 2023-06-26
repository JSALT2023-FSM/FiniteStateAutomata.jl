# SPDX-License-Identifier: CECILL-2.1

struct SparseSMMatrixCSC{K,Ti,TV<:AbstractVector{K},TI<:AbstractArray{Ti}} <: SparseSMMatrix{K}
    m::Int
    n::Int
    colptr::TI
    rowval::TI
    nzval::TV
end

#= AbstractArray API =#

Base.size(X::SparseSMMatrixCSC) = (X.m, X.n)

function Base.getindex(X::SparseSMMatrixCSC{K}, i::Int, j::Int) where K
    jj = X.colptr[j]:(X.colptr[j+1]-1)
    isempty(jj) && return zero(K)
    row = X.rowval[jj]
    val = X.nzval[jj]
    ii = findfirst(a -> a==i, row)
    (isnothing(ii)) ? zero(K) : val[ii]
end

function Base.getindex(X::SparseSMMatrixCSC{K}, ::Colon, j::Int) where K
    jj = getcolptr(X)[j]:(getcolptr(X)[j+1]-1)
    row = rowvals(X)[jj]
    val = nonzeros(X)[jj]
    SparseVector(size(X, 1), row, val)
end

Base.similar(X::SparseSMMatrixCSC) =
    SparseSMMatrixCSC(size(X)..., getcolptr(X), rowvals(X), similar(nonzeros(X)))
Base.similar(X::SparseSMMatrixCSC, T::Type) =
    SparseSMMatrixCSC(size(X)..., getcolptr(X), rowvals(X), similar(nonzeros(X), T))

#= Sparse API =#

getcolptr(X::SparseSMMatrixCSC) = getfield(X, :colptr)
rowvals(X::SparseSMMatrixCSC) = getfield(X, :rowval)
nonzeros(X::SparseSMMatrixCSC) = getfield(X, :nzval)

function _rangeptr(J, n)
    counts = fill!(similar(J, n+1), 0)
    counts[1] = 1
    for j in J
        counts[j+1] += 1
    end
    counts |> cumsum
end

#FIXME: Make sure for each column, the rows are sorted
function sparse(I::AbstractVector, J::AbstractVector, V::AbstractVector, m = maximum(I), n = maximum(J))
    p = sortperm(J)
    rowval = I[p]
    nzval = V[p]
    colptr = _rangeptr(J, n)
    SparseSMMatrixCSC(m, n, colptr, rowval, nzval)
end

sparse(I::AbstractVector, J::AbstractVector, val, m = maximum(I), n = maximum(J)) =
    sparse(I, J, repeat([val], length(I)), m ,n)

sparse(X::SparseSMMatrix) = sparse(findnz(X)..., X.m, X.n)

#FIXME: as it stands it will convert it to a non-CU matrix, ideally it should keep the proper CUDA type
#(without having to resort to something like cu(sparse(to_cpu(M))) )
function SparseSemimodules.sparse(M::AbstractMatrix{K}) where K
    I = Int64[]
    J = Int64[]
    V = K[]

	for i in 1:size(M, 2)
		v = M[:, i]
		nzind = findall(a->a!=zero(eltype(M)), v)
		nzval = v[nzind]
		append!(I, nzind)
		append!(J, fill(i, length(nzval)))
		append!(V, nzval)
	end

    if(M isa CuArray)
        return cu(SparseSemimodules.sparse(I, J, V, size(M)...))
    else
        return SparseSemimodules.sparse(I, J, V, size(M)...)
    end
end

function spdiagm(d::AbstractVector)
    N = length(d)
    colptr = similar(d, Int, N+1)
    colptr[1:end] = 1:N+1
    rowvals = similar(d, Int, N)
    rowvals[1:end] = 1:N
    SparseSMMatrixCSC(N, N, colptr, rowvals, d)
end

spzeros(::Type{K}, m::Integer, n::Integer) where K = sparse(Int[], Int[], zero(K), m, n)

nnz(X::SparseSMMatrixCSC) = length(X.nzval)

nnz(X::SparseSMMatrixCSC, i::Int) = X.colptr[i+1]-X.colptr[i]

function findnz(X::SparseSMMatrixCSC)
    J = similar(X.rowval, 0)
    for c in 1:size(X, 2)
        row = X[:, c]
        J = vcat(J, repeat([c], nnz(row)))
    end
    X.rowval, J, X.nzval
end

#= Concatenations =#
#TODO: Add the cases where xs is a mixed list of vector (normal or transpose) and matrices
#FIXME: The way we handle the different possibilities (mixed with csr, transpose) is not very good

Base.hcat(xs::Vararg{SparseSMVector}) = _hcat_spvector(xs)
Base.reduce(::typeof(hcat), xs::AbstractVector{<:SparseSMVector}) = _hcat_spvector(xs)

function _hcat_spvector(xs::Union{AbstractVector{<:SparseSMVector}, NTuple{N,<:SparseSMVector}}) where N
    m = length(xs[1])
    n = length(xs)
    rowval = reduce(vcat, [x.nzind for x in xs])
    nzval = reduce(vcat, [x.nzval for x in xs])

    counts = similar(rowval, n + 1)
    counts[1] = 1
    counts[2:end] = [nnz(x) for x in xs]
    colptr = cumsum(counts)

    SparseSMMatrixCSC(m, n, colptr, rowval, nzval)
end

# Horizontal

Base.hcat(Xs::Vararg{SparseSMMatrixCSC}) = _hcat_spmatrix(Xs)

Base.hcat(Xs::Vararg{<:SparseSMMatrix}) =
    _hcat_spmatrix(map(X -> (typeof(X)<:SparseSMMatrixCSR) ? sparse(X) : X, Xs))

Base.reduce(::typeof(hcat), Xs::AbstractVector{<:SparseSMMatrixCSC}) = _hcat_spmatrix(Xs)

function _hcat_spmatrix(Xs::Union{AbstractVector{<:SparseSMMatrixCSC}, NTuple{N,<:SparseSMMatrixCSC}}) where N
    dim_l = [size(X, 1) for X in Xs]
    all(dim_l .== size(Xs[1], 1)) || throw(DimensionMismatch("number of rows must be identical : $(dim_l)"))
    new_dims = (size(Xs[1])[1], sum([size(M)[2] for M in Xs]))

    I = vcat([rowvals(X) for X in Xs]...)
    colptr = reduce([getcolptr(X) for X in Xs], init=fill!(similar(getcolptr(Xs[1]), 1), 1)) do acc, c
        vcat(acc[1:end-1], c .+ (acc[end]-1))
    end
    V = vcat([nonzeros(X) for X in Xs]...)

    SparseSMMatrixCSC(new_dims..., colptr, I, V)
end

# Vertical

Base.vcat(Xs::Vararg{SparseSMMatrixCSC}) = _vcat_spmatrix(Xs)

Base.vcat(Xs::Vararg{<:SparseSMMatrix}) =
    _vcat_spmatrix(map(X -> (typeof(X)<:SparseSMMatrixCSR) ? sparse(X) : X, Xs))

function _vcat_spmatrix(Xs::Union{AbstractVector{<:SparseSMMatrixCSC}, NTuple{N,<:SparseSMMatrixCSC}}) where N
    dim_l = [size(X, 2) for X in Xs]
    all(dim_l .== size(Xs[1], 2)) || throw(DimensionMismatch("number of columns must be identical : $(dim_l)"))

    m_l = [size(M)[1] for M in Xs]
    tuple_l = [findnz(X) for X in Xs]

    new_dims = (sum(m_l), size(Xs[1])[2])
    I = []

    push_n = 0
    for i in eachindex(tuple_l) #[tᵢ[1] for tᵢ in tuple_l]
        n_row = push_n .+ tuple_l[i][1]
        push_n += m_l[i]
        push!(I, n_row)
    end

    I = collect(Iterators.flatten(I))
    J = collect(Iterators.flatten([tᵢ[2] for tᵢ in tuple_l]))
    V = collect(Iterators.flatten([tᵢ[3] for tᵢ in tuple_l]))

    return sparse(I, J, V, new_dims...)
end

# Horizontal and vertical

Base.hvcat(rows::Tuple{Vararg{Int64}}, Xs::Vararg{SparseSMMatrixCSC}) = _hvcat_spmatrix(rows, Xs)

Base.hvcat(rows::Tuple{Vararg{Int64}}, Xs::Vararg{<:SparseSMMatrix}) =
    _hvcat_spmatrix(rows, map(X -> (typeof(X)<:SparseSMMatrixCSR) ? sparse(X) : X, Xs))

function _hvcat_spmatrix(rows::Tuple{Vararg{Int64}}, Xs::Union{AbstractVector{<:SparseSMMatrixCSC}, NTuple{N,<:SparseSMMatrixCSC{K, Ti}}}) where {N, K, Ti}
    scan_l = fill!(Array{eltype(rows), 1}(undef, length(rows)), 0)
    for i in eachindex(rows) scan_l[i:end] .+= rows[i] end

    scan_l = [0;scan_l]

    buff_row = SparseSMMatrixCSC{K, Ti}[]

    for i in 2:length(scan_l)
        first = scan_l[i-1]+1
        last = scan_l[i]

        push!(buff_row, reduce([X for X in Xs[first+1:last]], init=Xs[first]) do acc, m
            hcat(acc, m)
        end)
    end

    return reduce(buff_row[2:end], init=buff_row[1]) do acc, r
        vcat(acc, r)
    end
end

#TODO: Generalize to both kind of matrices

function blockdiag(X::SparseSMMatrixCSC{Tv, Ti, TV, TI}...) where {Tv, Ti, TV, TI}

    num = length(X)
    mX = Int[ size(x, 1) for x in X ]
    nX = Int[ size(x, 2) for x in X ]
    m = sum(mX)
    n = sum(nX)

    colptr = TI(undef, n+1)
    nnzX = Int[ nnz(x) for x in X ]
    nnz_res = sum(nnzX)
    rowval = TI(undef, nnz_res)
    nzval = TV(undef, nnz_res)

    nnz_sofar = 0
    nX_sofar = 0
    mX_sofar = 0
    for i = 1 : num
        colptr[(1 : nX[i] + 1) .+ nX_sofar] = getcolptr(X[i]) .+ nnz_sofar
        rowval[(1 : nnzX[i]) .+ nnz_sofar] = rowvals(X[i]) .+ mX_sofar
        nzval[(1 : nnzX[i]) .+ nnz_sofar] = nonzeros(X[i])
        nnz_sofar += nnzX[i]
        nX_sofar += nX[i]
        mX_sofar += mX[i]
    end
    colptr[n+1] = nnz_sofar + 1

    SparseSMMatrixCSC(m, n, colptr, rowval, nzval)
end

#= Adapt =#

function Adapt.adapt_structure(to::Type{<:AbstractArray}, A::SparseSMMatrixCSC{Tv, Ti}) where {Tv, Ti}
	SparseSMMatrixCSC(size(A)..., adapt(findbasetype(to){Ti}, getcolptr(A)), adapt(findbasetype(to){Ti}, rowvals(A)), adapt(to, nonzeros(A)))
end

CUDA.cu(A::SparseSMMatrixCSC) = adapt(CuArray, A)
