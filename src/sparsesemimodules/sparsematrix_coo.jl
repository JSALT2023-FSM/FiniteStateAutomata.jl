import SparseArrays.sparse

abstract type AbstractSparseMatrixCOO{Tv, Ti <: Integer} <: AbstractSparseMatrix{Tv, Ti} end

abstract type AbstractSparseMatrix4KronCOO{S<:Semiring, Ti <: Integer} <: AbstractSparseMatrixCOO{S, Ti} end

macro lencheck(l, vars...)
    exprs = Expr[]
    for var in vars
      varname = string(var)
      push!(exprs, :(
        if length($(esc(var))) != $(esc(l))
          throw(DimensionError($varname, $(esc(l)), length($(esc(var)))))
    end
      ))
    end
    Expr(:block, exprs...)
  end

mutable struct SparseMatrixCOO{Tv, Ti <: Integer} <: AbstractSparseMatrixCOO{Tv, Ti}
    m::Int
    n::Int
    rows::Vector{Ti}
    cols::Vector{Ti}
    vals::Vector{Tv}

    function SparseMatrixCOO{Tv, Ti}(
    m::Integer,
    n::Integer,
    rows::Vector{Ti},
    cols::Vector{Ti},
    vals::Vector{Tv},
    ) where {Tv, Ti <: Integer}
    @noinline throwsz(str, lbl, k) =
        throw(ArgumentError("number of $str ($lbl) must be ≥ 0, got $k"))
    m < 0 && throwsz("rows", 'm', m)
    n < 0 && throwsz("columns", 'n', n)
    nnz = length(vals)
    # @lencheck nnz rows cols
    
    if !(nnz == length(rows) == length(cols))
      print(nnz,  length(rows), length(cols) )
      throw(AssertionError("lengths do not match"))
    end
    new(Int(m), Int(n), rows, cols, vals)
    end
end

function SparseMatrixCOO(m::Integer, n::Integer, rows::Vector, cols::Vector, vals::Vector)
    Tv = eltype(vals)
    Ti = promote_type(eltype(rows), eltype(cols))
    SparseArrays.sparse_check_Ti(m, n, Ti)
    nnz = length(vals)
    @lencheck nnz rows cols
    # silently shorten rowval and nzval to usable index positions.
    maxlen = abs(widemul(m, n))
    isbitstype(Ti) && (maxlen = min(maxlen, typemax(Ti) - 1))
    length(rows) > maxlen && resize!(rows, maxlen)
    length(cols) > maxlen && resize!(cols, maxlen)
    length(vals) > maxlen && resize!(vals, maxlen)
    SparseMatrixCOO{Tv, Ti}(m, n, rows, cols, vals)
  end

Base.size(A::SparseMatrixCOO) = (getfield(A, :m), getfield(A, :n))
        
# using matrix multiplication with csc, returns coo
# missing matmul in csr from semimodules 
function Base.:*(A::SparseMatrixCOO, B::SparseMatrixCOO)
	tocoo(sparse(A.rows,A.cols,A.vals, A.m, A.n)*sparse(B.rows,B.cols,B.vals, B.m, B.n))
end

function Base.:+(A::SparseMatrixCOO, B::SparseMatrixCOO)
	tocoo(sparse(A.rows,A.cols,A.vals, A.m, A.n)+sparse(B.rows,B.cols,B.vals, B.m, B.n))
end

# perform the kronocker product for spetial case used in composition
function Base.:*(A::SparseMatrixCOO{<:AbstractSparseMatrix4KronCOO,<:Int}, B::SparseMatrixCOO{<:AbstractSparseMatrix4KronCOO,<:Int})
	kron(A,B)
end

# TODO: complete
# function Base.show(io::IO, ::MIME"text/plain", S::AbstractSparseMatrixCOO)
#   	xnnz = nnz(S)
#   m, n = size(S)
#   print(io, m, "×", n, " ", typeof(S), " with ", xnnz, " stored ", xnnz == 1 ? "entry" : "entries")
#   if xnnz != 0
#     print(io, ":\n")
#     show(IOContext(io, :typeinfo => eltype(S)), S)
#   end
# end

function Base.print_array(io::IO, S::AbstractSparseMatrixCOO)
    if max(size(S)...) < 100
        D = todense(S)
        Base.print_matrix(io, D)
    else
        _show_with_braille_patterns(io, S)
    end
end

nnz(A::SparseMatrixCOO) = length(A.vals)
# nnz(A::Nothing) = 0

function Base.getindex(A::AbstractSparseMatrixCOO{Tv, Ti}, i0::Integer, i1::Integer) where {Tv, Ti}
    m, n = size(A)
    (1 ≤ i0 ≤ m && 1 ≤ i1 ≤ n) || throw(BoundsError())
    for k = 1:nnz(A)
        i, j, v = A.rows[k], A.cols[k], A.vals[k]
        if i == i0 && j == i1
            return v
        end
    end
end

hasitem(A::AbstractSparseMatrixCOO{Tv, Ti}, i0::Integer, i1::Integer) where {Tv, Ti} = typeof(A[i0,i1])!=Nothing