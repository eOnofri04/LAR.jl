LarA = Lar.Arrangement

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

function spatial_arrangement(model::Lar.Model)::Lar.Model

end

function pairwise_decomposition(model::Lar.Model, mt::Bool; err=1e-7)::Lar.Model
    sp_idx = LarA.spaceIndex(model, 2)
    de_models = [
        Lar.uniteModels!(
            de_model,
            face_decomposition(model, face_idx, sp_idx[face_idx])
        )
        for face_idx = 1 : size(model, 2, 1)
    ]
    return Lar.uniteModels(de_models, true, err)
end
