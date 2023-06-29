@testset "Vector FST" begin
    K = ProbSemiring{Float32}
    symtab = Dict(
        "a" => 1,
        "b" => 2,
        "c" => 3,
    )
    _arcs = [
        [(2, symtab["a"], symtab["a"], one(K))],
        [(2, symtab["b"], symtab["b"], one(K)), (3, symtab["c"], symtab["c"], one(K))],
        [(3, symtab["c"], symtab["c"], one(K))],
    ]
    finalweights = [zero(K), zero(K), one(K)]
    vfst = VectorFST(
        _arcs,
        1,
        finalweights
    )

    @test isinit(vfst, 1)
    @test isfinal(vfst, 3)
    @test !isfinal(vfst, 1) && !isfinal(vfst, 2)
    @test numstates(vfst) == 3
    @test numarcs(vfst, 2) == 2
    @test numarcs(vfst, 1) == 1

    q = addstate!(vfst)
    a_1qb = addarc!(vfst, 1, (q, symtab["b"], symtab["b"], one(K)))
    a_q3c = addarc!(vfst, q, (3, symtab["c"], symtab["c"], one(K)))

    @test numstates(vfst) == 4
    @test numarcs(vfst, 1) == 2
    @test numarcs(vfst, 2) == 2
    @test numarcs(vfst, 3) == 1

    @test length(finalstates(vfst)) == 1
    setfinalstate!(vfst, q)
    @test Set(finalstates(vfst)) == Set([3, q])

    deletearc!(vfst, 1, a_1qb)
    deletearcs!(vfst, q)
    @test numarcs(vfst, 1) == 1
    @test numarcs(vfst, q) == 0

    vfst = VectorFST(
        [
            [(2, 1, 1, K(.5)), (3, 2, 2, K(1.5))],
            [(3, 3, 3, K(2.5))],
            Arc{K}[]
        ],
        1,
        K[zero(K), zero(K), K(3.5)]
    )
    deletestate!(vfst, 2)
    @test numstates(vfst) == 2
    @test finalstates(vfst) == [2]
    @test sum([arcs(vfst, q) |> collect |> length for q in states(vfst)]) == 1


    @test_throws Exception deletestate!(vfst, 1) # 1 is init
    @test_throws Exception VectorFST(_arcs, 1, K[])
end
