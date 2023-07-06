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


@testset "Dense Composition" begin
    S = TropicalSemiring{Float32}

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
    A = convert(TensorFST{S, Array{S,4}}, A)
    B = convert(TensorFST{S, Array{S,4}}, B)
    C = dense_composition_sfo(A, B)

    coo = findall(x->x!=zero(S),M(C))
    print(coo)
    # @test coo == [CartesianIndex(5, 12, 1, 2),
    #         CartesianIndex(6, 12, 1, 2),
    #         CartesianIndex(2, 9, 2, 2),
    #         CartesianIndex(3, 9, 2, 2),
    #         CartesianIndex(5, 6, 3, 2),
    #         CartesianIndex(6, 6, 3, 2),
    #         CartesianIndex(1, 5, 1, 3),
    #         CartesianIndex(7, 11, 2, 3)]

    @test coo == [CartesianIndex(6, 12, 1, 2), 
    CartesianIndex(10, 12, 1, 2), 
    CartesianIndex(5, 11, 2, 2), 
    CartesianIndex(9, 11, 2, 2), 
    CartesianIndex(6, 10, 3, 2), 
    CartesianIndex(10, 10, 3, 2), 
    CartesianIndex(1, 6, 1, 3), 
    CartesianIndex(3, 8, 2, 3)]
end



@testset "Sparse States First Composition" begin
    S = TropicalSemiring{Float32}
    nsymbols = 3
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

    C = sparse_composition_sfo(A, B, nsymbols)

    truth = Vector{Tuple{Int64, Int64, Int64, TropicalSemiring{Float32}}}[
        [(5, 1, 3, S(0.4))],
        [(9, 2, 2, S(0.6))],
        [(9, 2, 2, S(0.8))],
        [],
        [(12, 1, 2, S(0.8)), (6, 3, 2, S(0.70000005))],
        [(12, 1, 2, S(1.0)), (6, 3, 2, S(0.90000004))],
        [(11, 2, 3, S(0.8))],
        [],
        [],
        [],
        [],
        []
        ]

    @test C.arcs == truth
end

@testset "Sparse Labels First Composition" begin
    S = TropicalSemiring{Float32}
    nsymbols = 3
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

    C = sparse_composition_lfo(A, B, nsymbols)
    print(C.arcs)
    
    truth = Vector{Tuple{Int64, Int64, Int64, TropicalSemiring{Float32}}}[
        [(5, 1, 3, S(0.4))],
        [(9, 2, 2, S(0.6))],
        [(9, 2, 2, S(0.8))],
        [],
        [(12, 1, 2, S(0.8)), (6, 3, 2, S(0.70000005))],
        [(12, 1, 2, S(1.0)), (6, 3, 2, S(0.90000004))],
        [(11, 2, 3, S(0.8))],
        [],
        [],
        [],
        [],
        []
        ]

    @test C.arcs == truth
    
end

