using SparseArrays
include("../src/Lar.jl")
LarA = Lar.Arrangement

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
σ = convert(Lar.Chain, SparseArrays.sparse([
    1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0
]))

model, sigma, edge_map = LarA.planar_arrangement_1(triforce, σ, true, false)
cleanedM, cleanedEM = LarA.cleandecomposition(model, sigma, edge_map)
bicon_comps = Lar.Arrangement.biconnected_components(cleanedM.T[1])
arrangedM = LarA.planar_arrangement_2(cleanedM, bicon_comps)
