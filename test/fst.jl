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
    @test numarcs(vfst, 1) == 1
    @test numarcs(vfst, q) == 1
    deletearcs!(vfst, q)
    @test numarcs(vfst, q) == 0

    setinitstate!(vfst, 2)
    @test isinit(vfst, 2)

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

    vfst = VectorFST(
        [
            [(2, 1, 1, K(.5)), (3, 2, 2, K(1.5))],
            [(3, 3, 3, K(2.5))],
            Tuple{Int,Int,Int,K}[]
        ],
        1,
        K[zero(K), zero(K), K(3.5)]
    )
    dfst = FiniteStateAutomata.densefst(vfst)
    @test size(M(dfst)) == (3, 3, 3, 3)
    dfst = convert(TensorFST{K, Array{K, 4}}, dfst)
    @test size(M(dfst)) == (3, 3, 3, 3)
end

@testset "Tensor FST" begin
    K = ProbSemiring{Float32}
    M = zeros(K, 3, 3, 2, 2)
    M[1,2,1,1] = one(K)
    M[2,3,2,2] = one(K)
    M[2,2,1,1] = one(K)
    α = [one(K), zeros(K, 2)...]
    ω = [zeros(K, 2)..., one(K)]
    dfst = TensorFST(M, α, ω)

    @test initstate(dfst) == 1
    @test Set(states(dfst)) == Set(1:3)
    @test numstates(dfst) == 3
    @test numarcs(dfst, 1) == 1 
    @test numarcs(dfst, 2) == 2
    @test arcs(dfst, 1) == [(2, 1, 1, one(K))]
end
