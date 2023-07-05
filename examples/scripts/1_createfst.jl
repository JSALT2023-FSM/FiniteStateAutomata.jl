
using ImageInTerminal
using Sixel


using FiniteStateAutomata

# Declare the semiring we want to use.
S = LogSemiring{Float32,1}

symbols = open(loadsymbols, "../data/symboltables/emoticons.syms")

# Compile a text-formatted FST
fst1 = compile(
    """
    1 2 1 1 -1
    2 2 2 2 0
    2 3 3 3 1
    1 3
    2 3.4
    """;
    semiring=S,
)

println("1. FST compiled from a string")
draw(fst1; isymbols=symbols, osymbols=symbols) |> dot(:png) |> rawdata -> display(MIME("image/png"), rawdata)

# Saving the dot format into a file.
open("out.dot", "w") do f
    write(f, draw(fst1))
end

# On the notebook
# draw(fst1) |> dot(:svg) |> HTML

println("2. FST compiled from a string compatible with OpenFST")
fst2 = compile(
    """
    0 1 1 1 -1
    1 1 2 2 0
    1 2 3 3 1
    0 2
    1 3.4
    """;
    semiring=S,
    openfst_compat=true
)

println("2. FST compiled from a string with 0-based index")
draw(fst1) |> dot(:png) |> rawdata -> display(MIME("image/png"), rawdata)

path = "../data/fsts/simple.txt"
fst3 = open(path, "r") do f
    compile(f; semiring=S)
end

println("3. FST compiled from a $path")
draw(fst3) |> dot(:png) |> rawdata -> display(MIME("image/png"), rawdata)

