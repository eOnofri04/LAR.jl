module Arrangement
    using IntervalTrees
    using NearestNeighbors
    using Triangle
	using SparseArrays
	using LinearAlgebra
	using Distributed
	using DataStructures
	using ..Lar;
	using LinearAlgebraicRepresentation

	LA = LinearAlgebra

    include("./minimal_cycles.jl")
	include("./spatial_indices.jl")
    # include("./dimension_travel.jl")
    include("./planar_arrangement.jl")
    include("./spatial_arrangement.jl")

    # export planar_arrangement, spatial_arrangement

end
