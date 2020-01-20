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


#-------------------------------------------------------------------------------
#   SPATIAL_INDEXING METHODS
#
#
#-------------------------------------------------------------------------------

"""
	buildBBIntervalTree(
		bboxes::Array{Array{Float64,2},1},
		dim::Int
	)::IntervalTrees.IntervalMap{Float64, Array}

Build the Interval Tree of a set of bounding boxes along the `dim` coordinate.
"""
function buildBBIntervalTree(
		bboxes::Array{Array{Float64,2},1},
		dim::Int
	)::IntervalTrees.IntervalMap{Float64, Array}

	# Build the dictionary related to the `dim`-coordinate
	boxDict = OrderedDict{Array{Float64,1}, Array{Int,1}}()
	for (h, box) in enumerate(bboxes)
		key = box[dim, :]
		if haskey(boxDict, key) == false
			boxDict[key] = [h]
		else
			push!(boxDict[key], h)
		end
	end

	# Build the related Interval Tree
	iT = IntervalTrees.IntervalMap{Float64, Array}()
	for (key, boxSet) in boxDict
		iT[tuple(key...)] = boxSet
	end

	return iT
end

"""
	boxcovering(
		bboxes::Array{Array{Float64,2},1},
		dim::Int
	)::Array{Array{Int,1},1}

Evaluates possible `bboxes` intersections over the `dim`-th coordinate

The function evaluates the pairwise intersection between the `dim` coordinate
of each `bboxes` couple. The operation is performed via Interval Tree.
"""
function boxcovering(
		bboxes::Array{Array{Float64,2},1},
		dim::Int
	)::Array{Array{Int,1},1}

	# Build the related Interval Tree
	tree = LarA.buildBBIntervalTree(bboxes, dim)

	# Build the covers
	covers = [[] for k = 1 : length(bboxes)]
	for (i, bbox) in enumerate(bboxes)
		extent = bboxes[i][dim, :]
		iterator = IntervalTrees.intersect(tree, tuple(extent...))
		for x in iterator
			append!(covers[i], x.value)
		end
	end

	return covers
end

#-------------------------------------------------------------------------------
#   SPATIAL_INDEXING
#-------------------------------------------------------------------------------

"""
	spaceIndex(model::Lar.Model)::Array{Array{Int,1},1}

Generation of *space indexes* for all ``(d-1)``-dim cell members of `model`.

*Spatial index* made by ``d`` *interval-trees* on
bounding boxes of ``sigma in S_{d−1}``. Spatial queries solved by
intersection of ``d`` queries on IntervalTrees generated by
bounding-boxes of geometric objects (LAR cells).

The return value is an array of arrays of `int`s, indexing cells whose
containment boxes are intersecting the containment box of the first cell.
According to Hoffmann, Hopcroft, and Karasick (1989) the worst-case complexity of
Boolean ops on such complexes equates the total sum of such numbers.

# Examples 2D

```
julia> V = hcat([[0.,0],[1,0],[1,1],[0,1],[2,1]]...);

julia> EV = [[1,2],[2,3],[3,4],[4,1],[1,5]];

julia> Sigma = Lar.spaceindex((V,EV))
5-element Array{Array{Int64,1},1}:
 [4, 5, 2]
 [1, 3, 5]
 [4, 5, 2]
 [1, 3, 5]
 [4, 1, 3, 2]
```

From `model2d` value, available in `?input_collection` docstring:

```julia
julia> Sigma =  spaceindex(model2d);
```

# Example 3D

```julia
model = model3d
Sigma =  spaceindex(model3d);
Sigma
```
"""
function spaceIndex(model::Lar.Model)::Array{Array{Int,1},1}

	bboxes = Lar.getModelBoundingBox(model, 1)
	# Build Box Coverings for each dimension and intersecting them
	covers = LarA.boxcovering(bboxes, 1)
	for d = 2 : length(model)
		dcovers = LarA.boxcovering(bboxes, d)
		covers  = [intersect(pair...) for pair in zip(covers, dcovers)]
	end

	# Remove each cell from its cover
	for k = 1 : length(covers)
		covers[k] = setdiff(covers[k],[k])
	end

	return covers
end
