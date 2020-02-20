using SparseArrays
using LinearAlgebraicRepresentation
LAR = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

V1, (VV,EV,FV,CV) = LAR.simplex(3, true)
tetra1 = V1, EV, FV, CV
V2 = [
  0.0 0.5 0.1 0.1
  0.0 0.1 0.5 0.1
  0.0 0.1 0.1 0.5
]
tetra2 = V2, EV, FV, CV
twotetra = LAR.Struct([ tetra, tetra2 ])
V,EV,FV,CV = LAR.struct2lar(twotetra)
GL.VIEW([ GL.GLGrid(V,FV, GL.Point4d(1,1,1,0.2)) ]);
m = Lar.Model(V);
Lar.addModelCells!(m, 1, convert(LAR.ChainOp, LAR.coboundary_0(EV::LAR.Cells)))
Lar.addModelCells!(m, 2, LAR.coboundary_1(V, FV::LAR.Cells, EV::LAR.Cells))
Lar.viewModel(m,1)


sp_idx = LarA.spaceIndex(m, 2)

de_models = [
	LarA.face_decomposition(m, face_idx, sp_idx[face_idx])
	for face_idx = 1 : size(m, 2, 1)
]

mfrag = Lar.mergeMultipleModels(de_models)

triangulated_faces = LAR.triangulate(convert(Lar.Points, mfrag.G'), mfrag.T[1:2])
FVs = convert(Array{Lar.Cells}, triangulated_faces);
GL.VIEW( GL.GLExplode(mfrag.G,FVs,expl,expl,expl,99,1) );

mArranged = Lar.deepcopy(mfrag)
Lar.addModelCells!(mArranged, 3, LAR.Arrangement.minimal_3cycles(
		convert(Lar.Points, mfrag.G'), mfrag.T[1], mfrag.T[2])
)




V,CVs,FVs,EVs = LAR.pols2tria(mArranged.G, mArranged.T[1], mArranged.T[2], mArranged.T[3])

GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
meshes = GL.GLExplode(V,CVs[1:end],8,4,6,99,0.5);
push!( meshes, GL.GLFrame)
GL.VIEW( meshes );

#=
cop_EV = LAR.coboundary_0(EV::LAR.Cells);
cop_EW = convert(LAR.ChainOp, cop_EV);
cop_FE = LAR.coboundary_1(V, FV::LAR.Cells, EV::LAR.Cells);
W = convert(LAR.Points, V');

V, copEV, copFE, copCF = LAR.Arrangement.spatial_arrangement( W::LAR.Points, cop_EW::LAR.ChainOp, cop_FE::LAR.ChainOp)

EV = LAR.cop2lar(copEV)
FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
FV = [collect(Set(cat(EV[e] for e in FE[f]))) for f=1:length(FE)]
FV = convert(LAR.Cells, FV)
W = convert(LAR.Points, V')
WW = [[k] for k=1:size(W,2)]

GL.VIEW(GL.numbering(.5)((W,[WW,EV]) ));

triangulated_faces = LAR.triangulate(V, [copEV, copFE])
FVs = convert(Array{LAR.Cells}, triangulated_faces)
V = convert(LAR.Points, V')
GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99));

EVs = LAR.FV2EVs(copEV, copFE) # polygonal face fragments
GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
