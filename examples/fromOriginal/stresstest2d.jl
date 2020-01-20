using LinearAlgebraicRepresentation
LAR = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

function stresstest2d()
	function input_collection(data::Array)::LAR.LAR
		assembly = LAR.Struct(data)
		return LAR.struct2lar(assembly)
	end

	V,(_,EV,FV) = LAR.cuboidGrid([4,1],true);
	W,(_,EW,FW) = LAR.cuboidGrid([1,5],true);
	mycircle(r,n) = LAR.circle(r)([n])
	data2d1 = (V,EV)
	data2d2 = LAR.Struct([ LAR.t(2,2), LAR.r(pi/3), LAR.t(-1.5,-2.5), (W,EW) ])
	data2d3 = LAR.Struct([ LAR.t(2,2), mycircle(2.5,160) ])
	data2d4 = LAR.Struct([ LAR.t(3.5,3.5), mycircle(.25,160) ])
	data2d5a = LAR.Struct([ LAR.t(5,3.5), mycircle(.75,160) ])
	data2d5 = LAR.Struct([ LAR.t(5,3.5), mycircle(.5,160) ])
	data2d6 = LAR.Struct([ LAR.t(5,3.5), mycircle(.25,160) ])
	data = [ data2d1, data2d2, data2d3, data2d4, data2d5,  data2d5a, data2d6 ]
	model2d = input_collection( [ LAR.Struct(data), LAR.t(-pi/6,pi/3), LAR.Struct(data) ] )
	V,EV = model2d
	GL.VIEW([ GL.GLLines(V,EV, GL.COLORS[1]) ]);

	model = Lar.Model(V);
	Lar.addModelCells!(model, 1, EV)

	model = Lar.Arrangement.planar_arrangement(model)

	EVs = Lar.FV2EVs(model)
	FVs = Lar.triangulate2D(model)
	GL.VIEW(GL.GLExplode(model.G, FVs, 1.2,1.2,1.2,99,1));
	GL.VIEW(GL.GLExplode(model.G, EVs, 1.2,1.2,1.2,99,1));
	GL.VIEW(GL.GLExplode(model.G, EVs, 1.2,1.2,1.2, 1,1));
end
