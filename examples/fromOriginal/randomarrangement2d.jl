using LinearAlgebraicRepresentation
LAR = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

function randomarrangement2d(n::Int = 100)
	V, EV = LAR.randomcuboids(n, .4)

	model = Lar.Model(GL.normalize2(V,flag=true))
	Lar.addModelCells!(model, 1, EV)
	model = Lar.Arrangement.planar_arrangement(model)
	FVs = Lar.triangulate2D(model)
	EVs = Lar.FV2EVs(model) # polygonal face fragments

	GL.VIEW(GL.GLExplode(model.G, EVs, 1.2,1.2,1.2,1,1));
	GL.VIEW(GL.GLExplode(model.G, FVs, 1.2,1.2,1.2,99,1));
	GL.VIEW(GL.GLExplode(model.G, FVs, 1,1,1,99,1));
end
