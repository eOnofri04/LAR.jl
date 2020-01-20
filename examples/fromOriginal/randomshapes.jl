using LinearAlgebraicRepresentation
LAR = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

function randomshapes(n::Int=100)
	V, EV = LAR.randomcuboids(n, .4)
	V = GL.normalize2(V, flag=true)
	GL.VIEW([ GL.GLLines(V, EV, GL.COLORS[1]) ])

	W, EW = LAR.fraglines(1.5,1.5,1.5)((V,EV))
	GL.VIEW([ GL.GLLines(W, EW, GL.COLORS[1]) ])

	model = Lar.Model(V)
	Lar.addModelCells!(model, 1, EV)
	model = Lar.Arrangement.planar_arrangement(model)

	FVs = Lar.triangulate2D(model)
	GL.VIEW(GL.GLExplode(model.G, FVs, 1.2,1.2,1.2,99,1));
	GL.VIEW(GL.GLExplode(model.G, FVs, 1,1,1,99,1));
end
