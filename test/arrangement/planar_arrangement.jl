LarA = Lar.Arrangement

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 1 METHODS
#    frag_edge_channel
#    frag_edge
#    intersect_edges
#    merge_vertices!
#-------------------------------------------------------------------------------

@testset "Edge Fragmentation" begin
    model = Lar.Model(hcat([
        [2., 2.], [4., 2.], [3., 3.5], [1., 3.], [5., 3.], [1., 2.], [5., 2.]
    ]...))
    model.T[1] = SparseArrays.sparse(Array{Int8, 2}([
        [1 1 0 0 0 0 0] #1->1,2
        [0 1 1 0 0 0 0] #2->2,3
        [1 0 1 0 0 0 0] #3->1,3
        [0 0 0 1 1 0 0] #4->4,5
        [0 0 0 0 0 1 1] #5->6,7
    ]))
    bigPI = Lar.spaceindex((model.G, Lar.cop2lar(model.T[1])))

    @testset "intersect_edges" begin
        V = convert(Lar.Points, model.G')
        EV = model.T[1]

        inters1 = LarA.intersect_edges(V, EV[5, :], EV[1, :])
        inters2 = LarA.intersect_edges(V, EV[1, :], EV[4, :])
        inters3 = LarA.intersect_edges(V, EV[1, :], EV[2, :])

        @test inters1 == [([2. 2.], 1/4),([4. 2.], 3/4)]
        @test inters2 == []
        @test inters3 == [([4. 2.], 1)]
    end

    @testset "frag_edges" begin
        V1, EV1 = LarA.frag_edge(model, 1, bigPI)
        V4, EV4 = LarA.frag_edge(model, 4, bigPI)
        V5, EV5 = LarA.frag_edge(model, 5, bigPI)

        @test V1 == [2.0 2.0; 4.0 2.0; 2.0 2.0; 4.0 2.0]
        @test EV1 == SparseArrays.sparse(Array{Int8, 2}([1 1 0 0]))

        @test isapprox(V4, [1.0 3.0; 5.0 3.0; 2.666666666 3.0; 3.33333333 3.0])
        @test EV4 == SparseArrays.sparse(Array{Int8, 2}([
            [1 0 1 0]
            [0 0 1 1]
            [0 1 0 1]
        ]))

        @test V5 == [1.0 2.0; 5.0 2.0; 2.0 2.0; 2.0 2.0; 4.0 2.0; 4.0 2.0]
        @test EV5 == SparseArrays.sparse(Array{Int8, 2}([
            [1 0 0 1 0 0]
            [0 0 0 1 0 1]
            [0 1 0 0 0 1]
        ]))
    end

    @testset "frag_edge_channel" begin
                                                                                # TODO
    end
end

@testset "merge_vertices" begin
    n0 = 1e-12
    n1l = 1-1e-12
    n1u = 1+1e-12
    model = Lar.Model([
        n0 -n0  n0 -n0  n0 -n0  n0 -n0 n1u n1l n1u n1l n1u n1l n1u n1l
        n0  n0 -n0 -n0 n1u n1u n1l n1l n1u n1u n1l n1l  n0  n0 -n0 -n0
    ])
    model.T[1] = SparseArrays.sparse(Array{Int8, 2}([
        [1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1]
        [1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]
        [0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1]
    ]))
    Lar.Arrangement.merge_vertices!(model)

    @test model.G == [n0 n0; n0 n1u; n1u n1u; n1u n0]'
    @test Array(model.T[1]) == [1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1]
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT 1
#-------------------------------------------------------------------------------

@testset "Planar Arrangement 1" begin
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - CLEAN DECOMPOSITION
#-------------------------------------------------------------------------------

@testset "Clean Decomposition" begin
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - BICONNECTED COMPONENTS
#-------------------------------------------------------------------------------

@testset "Clean Decomposition" begin
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2 - COMPONENT GRAPH (TGW) METHODS
#    get_external_cycle
#    pre_containment_test
#    prune_containment_graph
#    transitive_reduction!
#-------------------------------------------------------------------------------

@testset "TGW Methods" begin
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2 METHODS
#    componentgraph (aka TGW)
#    cell_merging
#-------------------------------------------------------------------------------

@testset "Component Graph" begin
end

@testset "Cell Merging" begin
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2
#-------------------------------------------------------------------------------

@testset "Planar Arrangement 2" begin
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE
#-------------------------------------------------------------------------------

@testset "Planar Arrangement" begin
end

#-------------------------------------------------------------------------------
#   COMPLETE PLANAR ARRANGEMENT
#-------------------------------------------------------------------------------

@testset "Complete Pipeline" begin

    #=== TRIFORCE MODEL DEFINITION ===#
    triforce = Lar.Model([
        0.0  4.0  8.0  4.0  2.0  6.0  4.0 -4.0 12.0  4.0  3.0  5.0
        0.0  6.0  0.0  0.0  3.0  3.0 -6.0  6.0  6.0  1.0  2.5  2.5
    ])
    triforce.T[1] = SparseArrays.sparse(Array{Int8, 2}([
        [1 1 0 0 0 0 0 0 0 0 0 0]
        [0 1 1 0 0 0 0 0 0 0 0 0]
        [1 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 1 0 0 0 0 0 0 0]
        [0 0 0 0 1 1 0 0 0 0 0 0]
        [0 0 0 1 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 1 0 0 0 0]
        [0 0 0 0 0 0 0 1 1 0 0 0]
        [0 0 0 0 0 0 1 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 1 0]
        [0 0 0 0 0 0 0 0 0 0 1 1]
        [0 0 0 0 0 0 0 0 0 1 0 1]
    ]))
    sigma = convert(Lar.Chain, SparseArrays.sparse([
        1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0
    ]))

    #=== PLANAR ARRANGEMENT 1 ===#
    arranged_triforce = Lar.Model([
        0.0  4.0  2.0  8.0  6.0  4.0  4.0 -4.0 12.0  4.0  3.0  5.0
        0.0  6.0  3.0  0.0  3.0  0.0 -6.0  6.0  6.0  1.0  2.5  2.5
    ])
    arranged_triforce.T[1] = SparseArrays.sparse(Int8[
        1 0 1 0 0 0 0 0 0 0 0 0
        0 1 1 0 0 0 0 0 0 0 0 0
        0 1 0 0 1 0 0 0 0 0 0 0
        0 0 0 1 1 0 0 0 0 0 0 0
        1 0 0 0 0 1 0 0 0 0 0 0
        0 0 0 1 0 1 0 0 0 0 0 0
        0 0 1 0 0 1 0 0 0 0 0 0
        0 0 1 0 1 0 0 0 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0
        1 0 0 0 0 0 1 0 0 0 0 0
        1 0 0 0 0 0 0 1 0 0 0 0
        0 1 0 0 0 0 0 1 0 0 0 0
        0 1 0 0 0 0 0 0 1 0 0 0
        0 0 0 1 0 0 1 0 0 0 0 0
        0 0 0 1 0 0 0 0 1 0 0 0
        0 0 0 0 0 0 0 0 0 1 1 0
        0 0 0 0 0 0 0 0 0 0 1 1
        0 0 0 0 0 0 0 0 0 1 0 1
    ])
    arranged_sigma = collect(1:6)
    arranged_edge_map = Array{Int64,1}[
        [ 1,  2], [ 3,  4], [ 5,  6],
        [     7], [     8], [     9],
        [10, 11], [12, 13], [14, 15],
        [    16], [    17], [    18]
    ]

    a, b, c = LarA.planar_arrangement_1(triforce, σ, true, false)
    @test arranged_triforce == a
    @test arranged_sigma    == b
    @test arranged_edge_map == c

    #=== CLEAN DECOMPOSITION ===#
    cleaned_triforce = Lar.Model([
        0.0 4.0 2.0 8.0 6.0 4.0 4.0 3.0 5.0;
        0.0 6.0 3.0 0.0 3.0 0.0 1.0 2.5 2.5
    ])
    cleaned_triforce.T[1] = SparseArrays.sparse(Int8[
        1 0 1 0 0 0 0 0 0
        0 1 1 0 0 0 0 0 0
        0 1 0 0 1 0 0 0 0
        0 0 0 1 1 0 0 0 0
        1 0 0 0 0 1 0 0 0
        0 0 0 1 0 1 0 0 0
        0 0 1 0 0 1 0 0 0
        0 0 1 0 1 0 0 0 0
        0 0 0 0 1 1 0 0 0
        0 0 0 0 0 0 1 1 0
        0 0 0 0 0 0 0 1 1
        0 0 0 0 0 0 1 0 1
    ])
    cleaned_edge_map = Array{Int64,1}[
        [1, 2], [3, 4], [5, 6],
        [   7], [   8], [   9],
        [    ], [    ], [    ],
        [  10], [  11], [  12]
    ]

    a, b = LarA.cleandecomposition(a, b, c)
    @test cleaned_triforce == a
    @test cleaned_edge_map == b

    #=== BICONNECTED COMPONENTS ===#
    bicon_comps = Array{Int64,1}[[8, 9, 7, 5, 6, 4, 3, 2, 1], [12, 11, 10]]
    c = Lar.Arrangement.biconnected_components(a.T[1])
    @test bicon_comps == c

    #=== PLANAR ARRANGEMENT 2 ===#

    rearranged_triforce = Lar.Model([
        0.0 4.0 2.0 8.0 6.0 4.0 4.0 3.0 5.0
        0.0 6.0 3.0 0.0 3.0 0.0 1.0 2.5 2.5
    ])
    rearranged_triforce.T[1] = SparseArrays.sparse(Int8[
        -1  0  1  0  0  0  0  0  0
         0 -1  1  0  0  0  0  0  0
         0 -1  0  0  1  0  0  0  0
         0  0  0 -1  1  0  0  0  0
        -1  0  0  0  0  1  0  0  0
         0  0  0 -1  0  1  0  0  0
         0  0 -1  0  0  1  0  0  0
         0  0 -1  0  1  0  0  0  0
         0  0  0  0 -1  1  0  0  0
         0  0  0  0  0  0 -1  1  0
         0  0  0  0  0  0  0 -1  1
         0  0  0  0  0  0 -1  0  1
    ])
    rearranged_triforce.T[2] = SparseArrays.sparse([
        0  0  0  0  0  0  1 -1 -1  1  1 -1
        -1 0  0  0  1  0 -1  0  0  0  0  0
        0  1 -1  0  0  0  0  1  0  0  0  0
        0  0  0  1  0 -1  0  0  1  0  0  0
        0  0  0  0  0  0  0  0  0 -1 -1  1
    ])
    a = Lar.Arrangement.planar_arrangement_2(a, c, b)
    @test rearranged_triforce == a

end
