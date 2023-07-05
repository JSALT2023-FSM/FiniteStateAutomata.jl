
using FiniteStateAutomata

S = LogSemiring{Float32,1}

fst = open("../examples/fsts/fst_ex1.txt") do f
    compile(f; semiring = S)
end

display(MIME("text/vnd.graphviz"), draw(fst))
display(MIME("text/plain"), draw(fst))
display(MIME("image/svg+xml"), draw(fst))

