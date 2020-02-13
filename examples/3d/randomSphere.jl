@using LinearAlgebraicRepresentation as LAR
@using LinearAlgebra as LA
@using ViewerGL as GL
include("../../src/Lar.jl")
LarA = Lar.Arrangement

function build_copFE(FV::Lar.Cells, EV::Lar.Cells)
	copFE = LAR.u_coboundary_1(FV, EV) # unsigned
	faceedges = [findnz(copFE[f,:])[1] for f=1:size(copFE,1)]

	f_edgepairs = Array{Array{Int64,1}}[]
	for f=1:size(copFE,1)
		edgepairs = Array{Int64,1}[]
		for v in FV[f]
			push!(edgepairs, [e for e in faceedges[f] if v in EV[e]])
		end
		push!(f_edgepairs, edgepairs)
	end
	copFE.=0
	for f=1:size(copFE,1)
		facepairs = sort(sort.(f_edgepairs[f]))
		# setting a reference
		copFE[f, facepairs[1][1]] = 1
		for (e1,e2) in facepairs
			@assert copFE[f, e1] != 0
			if copFE[f, e2] == 0
				copFE[f, e2] =
					-((EV[e1][1]==EV[e2][1])||(EV[e1][2]==EV[e2][2])) * copFE[f,e1]+
					((EV[e1][1]==EV[e2][2])||(EV[e1][2]==EV[e2][1])) * copFE[f,e1]
			end
		end
	end
	return copFE
end

function getGFVModel(V::Lar.Points, FV::Array{Array{Int,1}})::Lar.Model
	function getFaceEdges(fV::Array{Int,1})
		return [fV[1], fV[2]], [fV[1], fV[3]], [fV[2], fV[3]]
	end
	model = Lar.Model(V);
	EV = unique([(map(f -> getFaceEdges(f), FV)...)...])
	Lar.addModelCells!(model, 1, LAR.build_copEV(EV))
	#FE = [[3x+1, 3x+2, 3x+3] for x = 0 : length(FV)-1]
	Lar.addModelCells!(model, 2, build_copFE(FV, EV))
	return Lar.mergeModelVertices(model)
end

scaling = 1.2
no_sphere = 4
expl = 2.0

carry = Array{Tuple{Array{Float64,2},Array{Array{Int64,1},1}}}(undef, no_sphere)
for i = 1 : no_sphere
	carry[i] = LAR.struct2lar(LAR.Struct([
		LAR.t(rand(0.5:1e-2:5, 3)...),
		LAR.sphere(rand(1:1e-2:3) * scaling)()
	]))
end
V, FV = LAR.struct2lar(LAR.Struct(carry))
FV = sort(sort.(FV))
GL.VIEW([GL.GLGrid(V, FV)])

m = getGFVModel(V, FV)
Lar.viewModel(m, 1)
Lar.viewModel(m, 2)
GL.VIEW(GL.GLExplode(m.G,[[s] for s in FV],1.2,1.2,1.2,99,1));

sp_idx = LarA.spaceIndex(m, 2)

de_models = [
	LarA.face_decomposition(m, face_idx, sp_idx[face_idx])
	for face_idx = 1 : size(m, 2, 1)
]

m = deepcopy(de_models[1])
for i = 2 : length(de_models)
	Lar.uniteModels!(m, de_models[i])
end

triangulated_faces = LAR.triangulate2D(convert(Lar.Points, m.G'), m.T[1:2])
FVs = convert(Array{Lar.Cells}, triangulated_faces);
GL.VIEW( GL.GLExplode(m.G,FVs,expl,expl,expl,99,1) );

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
