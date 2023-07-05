
using ImageInTerminal
using Sixel
using FiniteStateAutomata

# Declare the semiring we want to use.
S = LogSemiring{Float32,1}

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
draw(fst1) |> dot(:png) |> rawdata -> display(MIME("image/png"), rawdata)

println("2. FST compiled from a string compatible with OpenFST")

path = "../data/fsts/simple.txt"
fst3 = open(path, "r") do f
    # openfst_compat=true read the file as a 0-based state indexing.
    compile(f; semiring=S)
end

println("3. FST compiled from a $path")
draw(fst3) |> dot(:png) |> rawdata -> display(MIME("image/png"), rawdata)

