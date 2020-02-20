using SparseArrays, Test
LarA = Lar.Arrangement
using LinearAlgebraicRepresentation
LAR = LinearAlgebraicRepresentation


@testset "Planar Arrangement - Big Random Data" begin
    

    @testset "Closure of 2-cells" begin
        for i = 1 : size(model2, 2, 1)
            edges  = findnz(model2.T[2][i, :])[1]
            points = abs.(model2.T[1][edges, :])
            @test sum(points) == 2 * length(edges)
        end
    end

    @testset "Well formed cells" begin
        @test sum(sum(abs.(model2.T[1]), dims = 1) .<= 1) == 0
        @test sum(sum(abs.(model2.T[1]), dims = 2) .!= 2) == 0
        @test sum(sum(abs.(model2.T[2]), dims = 1) .>= 3) == 0
        @test sum(sum(abs.(model2.T[2]), dims = 2) .<= 2) == 0
    end

    @testset "Euler" begin
        @test size(model2, 0, 2) - size(model2, 1, 1) + size(model2, 2, 1) == 1
    end

    @testset "1-Cell Compactness" begin
        #@test sum(abs.(model2.T[1])) == 2 * size(model, 1, 2)
    end

    @testset "2-Cell Compactness" begin
        #@test sum(abs.(model2.T[2])) == 2 * size(model2, 2, 2)
    end

    @testset "Biconnected Components" begin
        # for comp in bicon_comps
        #     sum(
        #         (sum(abs.(model.T[1][comp, :]), dims=1) .!= 0)
        #         + (sum(abs.(model.T[1][setdiff(1:size(model, 1, 1), comp), :]), dims=1) .!= 0)
        #         .!= 1
        #     ) == 0
        # end
    end
end
