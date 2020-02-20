LarA = Lar.Arrangement
LAR = LinearAlgebraicRepresentation

using IntervalTrees
using SparseArrays
using NearestNeighbors
using DataStructures

#-------------------------------------------------------------------------------
#	TSAS PIPELINE
#    input_collection
#    spatial_indexing
#    pairwise_decomposition
#    equivalence_congruence
#    connection_extraction
#    topological_gift_wrapping
#    boundaries_evaluation
#    holes_cotainments
#    transitive_reduction
#    cell_adjoining
#    cell_assembling
#    colleclt_output
#-------------------------------------------------------------------------------

function spatial_arrangement(
		model::Lar.Model,
		mt::Bool = false;
		err=1e-7
	)::Lar.Model
    model = Lar.Arrangement.pairwise_decomposition(model, mt, err = err);
    bicon_comps = LarA.biconnected_components(model.T[1]);
    ##arr_models = [purgeHoles(TGW(...)) for bicon in bicon_comps] # via spaceIndex(model, 3)
    Lar.mergeMultipleModels(arr_models, true, err)
end

function pairwise_decomposition(
		model::Lar.Model,
		mt::Bool = false;
		err=1e-7
	)::Lar.Model

	!mt || throw(ArgumentError("Multicore not Coded!"))
    sp_idx = LarA.spaceIndex(model, 2)
    de_models = [
        LarA.face_decomposition(model, face_idx, sp_idx[face_idx])
        for face_idx = 1 : size(model, 2, 1)
    ]
    return Lar.mergeMultipleModels(de_models, err=err)
end

function face_decomposition(#V, EV, FE, sp_idx, sigma)
		model::Lar.Model,
		face_idx::Int,
		sp_idx_face::Array{Int,1}
	)::Lar.Model

	function build_projection_matrix(vs::Lar.Points)::Matrix{Float64}
    	u1 = vs[:, 2] - vs[:, 1]
    	u2 = vs[:, 3] - vs[:, 1]
    	u3 = LAR.cross(u1, u2)
    	T = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    	T[1:3, 4] = - vs[:,1]
    	M = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    	M[1:3, 1:3] = [u1 u2 u3]
    	return M' * T
	end

	function face_intersection(													# Inserire il modello in depth
			V::Lar.Points,
			EV::Lar.ChainOp,
			face::Lar.Cell
		)::Lar.Model

	    vs = LAR.buildFV(EV, face)
	    retV = Lar.Points(undef, 2, 0)

	    visited_verts = []
	    for i in 1:length(vs)
	        o = V[:, vs[i]]
	        j = i < length(vs) ? i+1 : 1
	        d = V[:, vs[j]] - o

	        err = 10e-8
	        if !(-err < d[3] < err)

	            alpha = -o[3] / d[3]

	            if -err <= alpha <= 1+err
	                p = o + alpha*d

	                if -err < alpha < err || 1-err < alpha < 1+err
	                    if !(LAR.vin(p, visited_verts))
	                        push!(visited_verts, p)
	                        retV = [retV reshape(p, 3, 1)[1:2, 1]]
	                    end
	                else
	                    retV = [retV reshape(p, 3, 1)[1:2, 1]]
	                end
	            end
	        end

	    end

	    vnum = size(retV, 2)

	    if vnum == 1
	        vnum = 0
	        retV = Lar.Points(undef, 2, 0)
	    end
	    enum = (÷)(vnum, 2)														##?
	    retEV = spzeros(Int8, enum, vnum)

	    for i in 1:enum
	        retEV[i, 2*i-1:2*i] = [-1, 1]
	    end

	    retModel = Lar.Model(retV)
		Lar.addModelCells!(retModel, 1, retEV)
		return retModel
	end

	#== THE FUNCTION ==#

	Gface, Gfaceidx = Lar.getModelCellGeometry(model, 2, face_idx, true)
	M = build_projection_matrix(Gface)
	PG = (M * [model.G; ones(1, size(model, 0, 2))])[1:3, :]					# Tutti i vertici?

	# Generate face Model (in 2D)
	Pmodel = Lar.Model(PG[1:2, Gfaceidx])
	Lar.addModelCells!(Pmodel, 1, model.T[1][Lar.getModelLoCell(model, 2, face_idx), Gfaceidx])

	# Each other face adds new pieces
	for f in sp_idx_face
		println("Intersecting $face_idx with $f")
		newModel = face_intersection(PG, model.T[1], model.T[2][f, :])
		Lar.uniteModels!(Pmodel, newModel)
	end

	Pmodel = Lar.mergeModelVertices(Pmodel)
													#---v ? is needed sigma?
	model = Lar.Arrangement.planar_arrangement(Pmodel, sparsevec(ones(Int8, length(Gfaceidx))))
	@assert !isnothing(model) "UNEXPECTED ERROR: a face should be mapped to itself"

	n = size(model, 0, 2)
	retModel = Lar.Model((inv(M)*[model.G; zeros(1, n); ones(1, n)])[1:3, :])
	Lar.addModelCells!(retModel, 1, model.T[1])
	Lar.addModelCells!(retModel, 2, model.T[2])

	return retModel
end

#	vs_num = size(V, 1)

	# 2D transformation of sigma face
#    sigmavs = (abs.(FE[sigma:sigma,:])*abs.(EV))[1,:].nzind
#    sV = V[sigmavs, :]
#    sEV = EV[FE[sigma, :].nzind, sigmavs]
#    M = Lar.Arrangement.submanifold_mapping(sV)
#    tV = ([V ones(vs_num)]*M)[:, 1:3]  # folle convertire *tutti* i vertici
#    sV = tV[sigmavs, :]
    # sigma face intersection with faces in sp_idx[sigma]
#    for i in sp_idx[sigma]
#        tmpV, tmpEV = Lar.Arrangement.face_int(tV, EV, FE[i, :])
#
#        sV, sEV = Lar.skel_merge(sV, sEV, tmpV, tmpEV)
#    end

    # computation of 2D arrangement of sigma face
#    sV = sV[:, 1:2]
#    nV, nEV, nFE = planar_arrangement(sV, sEV, sparsevec(ones(Int8, length(sigmavs))))
#    if nV == nothing ## not possible !! ... (each original face maps to its decomposition)
#        return [], spzeros(Int8, 0,0), spzeros(Int8, 0,0)
#    end
#    nvsize = size(nV, 1)
#    nV = [nV zeros(nvsize) ones(nvsize)]*inv(M)[:, 1:3] ## ????
#    return nV, nEV, nFE

#=
function pairwise_decomposition(model::Lar.Model, mt::Bool; err=1e-7)::Lar.Model
    sp_idx = LarA.spaceIndex(model, 2)
    de_models = [
        Lar.uniteModels!(
            de_model,
            LarA.face_decomposition(model, face_idx, sp_idx[face_idx])
        )
        for face_idx = 1 : size(model, 2, 1)
    ]
    return Lar.uniteModels(de_models, true, err)
end
=#
