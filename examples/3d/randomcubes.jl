@using LinearAlgebraicRepresentation as LAR
@using LinearAlgebra as LA
@using ViewerGL as GL

include("../../src/Lar.jl")

store = [];
scaling = 1.20;
V,(VV,EV,FV,CV) = LAR.cuboid([0.25,0.25,0.25],true,[-0.25,-0.25,-0.25]);
mybox = (V,CV,FV,EV);
for k=1:3
	size = rand()*scaling
	scale = LAR.s(size,size,size)
	transl = LAR.t(rand(3)...)
	alpha = 2*pi*rand()
	rx = LAR.r(alpha,0,0); ry = LAR.r(0,alpha,0); rz = LAR.r(0,0,alpha)
	rot = rx * ry * rz
	str = LAR.Struct([ transl, scale, rot, mybox ])
	obj = LAR.struct2lar(str)
	vs = obj[1]
	diag = LinearAlgebra.norm(vs[:,8]-vs[:,1])
	if diag > 1/5
		push!(store, obj)
	end
end

str = LAR.Struct(store);
V,CV,FV,EV = LAR.struct2lar(str);

GL.VIEW([ GL.GLPol(V,CV, GL.COLORS[1]) ]);

cop_EV = LAR.coboundary_0(EV::Lar.Cells);
cop_FE = LAR.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);

model = Lar.Model(V);
Lar.addModelCells!(model, 1, cop_EV)
Lar.addModelCells!(model, 2, cop_FE)

#-------------------------------------------------------------------------------
#
# V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement(
# 	 W::Lar.Points, cop_EV::Lar.ChainOp, cop_FE::Lar.ChainOp);
#
#
# triangulated_faces = Lar.triangulate2D(V, [copEV, copFE]);
# FVs = convert(Array{Lar.Cells}, triangulated_faces);
# V = convert(Lar.Points, V');
# GL.VIEW( GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1) );
#
#
# EVs = Lar.FV2EVs(copEV, copFE); # polygonal face fragments
# GL.VIEW( GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1) );
