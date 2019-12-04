#-------------------------------------------------------------------------------
#   MODEL STRUCT DEFINITION
#-------------------------------------------------------------------------------

mutable struct Model
    G::Points
    T::Array{ChainOp, 1}

    function Model(V::Points)
        size(V, 1) > 0 ||
            throw(ArgumentError("At least one point is needed."))
        new(V, Array{ChainOp, 1}(undef, size(V, 1)))
    end

    function Model(V::Points, T::Array{ChainOp, 1})
        dim, npts = size(V)
        dim > 0 ||
            throw(ArgumentError("At least one point is needed."))
        length(T) == dim ||
            throw(ArgumentError("Topology is not coherent with Geometry."))
        for i = 1 : dim
            !isdefined(T, i) || size(T, 2) == npts ||
                throw(ArgumentError("Topology is not coherent with Geometry."))
        end
        new(V, T)
    end

    function Model()
        nothing
    end
end

#-------------------------------------------------------------------------------
#   BASIC PROPERTIES
#-------------------------------------------------------------------------------

length(m::Lar.Model)    = size(m.G, 1)
size(m::Lar.Model)		= size(m.G)
Base.copy(m::Lar.Model) = Lar.Model(m.G, m.T)

#==function +(m1::Lar.Model, m2::Lar.Model)
    if m1.dim > m2.dim
        V1 = m1.G
        mdiff = m1.dim - m2.dim
        V2 = [m2.G; zeros(mdiff, m2.n)]
    elseif m1.dim < m2.dim
        mdiff = m2.dim - m1.dim
        V1 = [m1.G; zeros(mdiff, m1.n)]
        V2 = m2.G
    else
        V1 = m1.G
        V2 = m2.G
    end
    Lar.Model([V1 V2])
end

function checkModel(m::Model)
    @assert size(m.G) == (m.dim, m.n);
    @assert length(m.T) == m.dim;
    for i = 1 : m.dim
        @assert !isdefined(m.T, i) || length(m.T[i]) == n;
    end
end==#

#-------------------------------------------------------------------------------
#   GEOMETRY MANIPULATION
#-------------------------------------------------------------------------------

function addModelVertex!(m::Lar.Model, v::Array{Float64, 1})::Nothing
    length(V) == length(v) ||
        throw(ArgumentError("Point dimension mismatch."))
    m.G = [m.G v];
end

function deleteModelVertex!(m::Lar.Model, v::Int)::Nothing
    deleteModelVertices!(m, [v])
end

function deleteModelVertices!(m::Lar.Model, vs::Array{Int,1})::Nothing
    m.T[1] = m.T[1][:, setdiff(collect(1 : size(m.G, 2)), vs)]
    return
end

#-------------------------------------------------------------------------------
#   TOPOLOGY MANIPULATION
#-------------------------------------------------------------------------------

function getModelEdgeVertices(m::Lar.Model, e::Int)
    vidxs  = m.T[1][e, :].nzind
    v1, v2 = map(i->m.G[:, vidxs[i]], [1,2])
    return (v1, v2)
end

function deleteModelEdge!(m::Lar.Model, e::Int)::Nothing
    deleteModelEdges!(m, [e])
end

function deleteModelEdges!(m::Lar.Model, ev::Array{Int,1})::Nothing
    EVs     = m.T[1]
    EVs     = EVs[  setdiff(collect(1 : EVs.m), ev),  :  ]
    tokeepV = [length(EVs[:,i].nzval) != 0  for i = 1 : EVs.n]
    m.G     = m.G[:, tokeepV]
    m.T[1]  = EVs[:, tokeepV]
    return
end
