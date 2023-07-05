S = ProbSemiring{Float32}

fsttxt_origin = """1 2 1 1
    2 3 2 2
    3 4 3 3
    4 5 4 4 0.5
    5"""

for S in [ProbSemiring{Float32}, LogSemiring{Float32, 1}, TropicalSemiring{Float32}]
    fst = compile(fsttxt_origin; semiring = S)

    @test fst isa VectorFST
    @test numstates(fst) == 5
    @test numarcs(fst, 1) == 1 && first(arcs(fst, 1)) == (2, 1, 1, one(S))
    @test numarcs(fst, 2) == 1 && first(arcs(fst, 2)) == (3, 2, 2, one(S))
    @test numarcs(fst, 3) == 1 && first(arcs(fst, 3)) == (4, 3, 3, one(S))
    @test numarcs(fst, 4) == 1 && first(arcs(fst, 4)) == (5, 4, 4, S(0.5))
    @test numarcs(fst, 5) == 0
    @test finalweight(fst, 5) == one(S)
    @test initstate(fst) == 1


    buf = IOBuffer()
    print(IOContext(buf, :compact => true), fst)
    @test String(take!(buf)) == fsttxt_origin
end
