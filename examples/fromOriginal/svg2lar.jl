using LinearAlgebraicRepresentation
LAR = LinearAlgebraicRepresentation
using ViewerGL,SparseArrays
GL = ViewerGL

function svgarrangement(filename)
	V, EV = LAR.svg2lar(filename)
	GL.VIEW([ GL.GLLines(V,EV) ])
	model = Lar.Model(V)
	Lar.addModelCells!(model, 1, EV)
	model = Lar.Arrangement.planar_arrangement(model)
	bicon_comps = Lar.Arrangement.biconnected_components(model.T[1])

	# visualization of component graphs
	EV = LAR.cop2lar(model.T[1])
	V, EVs = LAR.biconnectedComponent((model.G, EV::LAR.Cells)) # 2-connected components
	VV = [[k] for k = 1 : size(V,2)]
	meshes = [ GL.numbering(.1)((V, (VV,EVs[k])), GL.COLORS[k],1) for k=1:length(EVs) ]
	GL.VIEW( cat(meshes) );

	# final solid visualizations
	# FE = [SparseArrays.findnz(model.T[2][k,:])[1] for k=1:size(model.T[2],1)]
	# FV = [collect(Set(cat([EV[e] for e in FE[f]]))) for f=1:length(FE)]

	EVs = Lar.FV2EVs(model)
	FVs = Lar.triangulate2D(model)
	GL.VIEW(GL.GLExplode(model.G, FVs, 1.2, 1.2, 1.2, 99, 1));
	GL.VIEW(GL.GLExplode(model.G, FVs,   1,   1,   1, 99, 1));
	GL.VIEW(GL.GLExplode(model.G, EVs,   1,   1,   1, 99, 1));
end
