@testset "Vector FST" begin
    K = ProbSemiring{Float32}
    symtab = Dict(
        "a" => 1,
        "b" => 2,
        "c" => 3,
    )
    states = [
        [(2, symtab["a"], one(K))],
        [(2, symtab["b"], one(K)), (3, symtab["c"], one(K))],
        [(3, symtab["c"], one(K))],
    ]
    finalweights = [zero(K), zero(K), one(K)]
    vfst = VectorFST(
        states,
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
    a_1qb = addarcs!(vfst, 1, (q, symtab["b"], one(K)))
    a_q3c = addarcs!(vfst, q, (3, symtab["c"], one(K)))

    @test numstates(vfst) == 4
    @test numarcs(vfst, 1) == 2
    @test numarcs(vfst, 2) == 2
    @test numarcs(vfst, 3) == 1

    deletearc!(vfst, 1, a_1qb)
    deletearcs!(vfst, q)
    @test numarcs(vfst, 1) == 1
    @test numarcs(vfst, q) == 0

    deletestate!(vfst, q)
    @test numstates(vfst) == 3
end
