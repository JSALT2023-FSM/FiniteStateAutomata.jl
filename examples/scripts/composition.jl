using SparseArrays
using FiniteStateAutomata
using ImageInTerminal
using Sixel

symbols = Dict(1 => "a", 2 => "b", 3 => "c")

S = TropicalSemiring{Float64}

A = VectorFST(
    [
        Arc{S}[(2, 1, 2, S(0.1)), (3, 2, 1, S(0.2))],
        Arc{S}[(2, 3, 1, S(0.3)), (4, 1, 1, S(0.4))],
        Arc{S}[(4, 2, 2, S(0.5))],
        Arc{S}[]
    ],
    1,
    [zero(S),zero(S),zero(S),S(0.6)]
)
B = VectorFST(
    [
        Arc{S}[(2, 2, 3, S(0.3))],
        Arc{S}[(3, 1, 2, S(0.4))],
        Arc{S}[(3, 1, 2, S(0.6))],
    ],
    1,
    S[zero(S),zero(S),S(0.7)]
)

draw(A; isymbols=symbols, osymbols=symbols) |> dot(:png) |> rawdata -> display(MIME("image/png"), rawdata)

draw(B; isymbols=symbols, osymbols=symbols) |> dot(:png) |> rawdata -> display(MIME("image/png"), rawdata)

C = sparse_composition_sfo(A, B, length(symbols))

draw(C; isymbols=symbols, osymbols=symbols) |> dot(:png) |> rawdata -> display(MIME("image/png"), rawdata)
