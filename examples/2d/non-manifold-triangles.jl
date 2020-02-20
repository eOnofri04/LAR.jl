@using LinearAlgebraicRepresentation as LAR
@using ViewerGL as GL
using SparseArrays
include("../../src/Lar.jl")
LarA = Lar.Arrangement

model = Lar.Model([
    0.0 1.0 0.0 0.5 0.1
    0.0 0.0 1.0 0.1 0.5
])
Lar.addModelCells!(model, 1, SparseArrays.sparse(Array{Int8, 2}([
    -1  1  0  0  0
    -1  0  1  0  0
     0 -1  1  0  0
    -1  0  0  1  0
    -1  0  0  0  1
     0  0  0 -1  1
])))
Lar.addModelCells!(model, 2, SparseArrays.sparse(Array{Int8, 2}([
     1 -1  1  0  0  0
     0  0  0  1 -1  1
])))
Lar.view2DModel(model)
arrangedModel = LarA.planar_arrangement(model)
Lar.view2DModel(arrangedModel, 1)
