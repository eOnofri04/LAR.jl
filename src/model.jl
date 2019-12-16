import Base: +, length, size, ==

#-------------------------------------------------------------------------------
#   MODEL STRUCT DEFINITION
#-------------------------------------------------------------------------------

mutable struct Model
    G::Lar.Points
    T::Array{Lar.ChainOp, 1}
    topological::Bool

    function Model(V::Lar.Points, T::Array{Lar.ChainOp, 1}, topological::Bool)
        dim, npts = size(V)
        dim > 0 ||
            throw(ArgumentError("At least one point is needed."))
        length(T) == dim ||
            throw(ArgumentError("Topology is not coherent with Geometry."))
        if topological
            isempty(T[1]) || size(T[1], 2) == npts ||
                throw(ArgumentError("Topology not coherent with Geometry."))
            for i = 2 : dim
                isempty(T[i-1]) || isempty(T[i]) ||
                    size(T[i], 2) == size(T[i-1],1) ||
                    throw(ArgumentError("Topology not coherent with Topology."))
            end
        else
            for i = 1 : dim
                isempty(T[i]) || size(T[i], 2) == npts ||
                    throw(ArgumentError("Topology not coherent with Geometry."))
            end
        end
        new(V, T, topological)
    end
    function Model(V::Lar.Points, T::Array{Lar.ChainOp, 1})
        Lar.Model(V, T, true)
    end
    function Model(V::Lar.Points, topological::Bool)
        nmat = !topological * size(V, 2);
        T = convert(Array{Lar.ChainOp,1},
            [SparseArrays.spzeros(Int8, 0, nmat) for i = 1 : size(V, 1)]
        )
        if topological
            T[1] = convert(Lar.ChainOp, SparseArrays.spzeros(Int8, 0,size(V,2)))
        end
        Lar.Model(V, T, topological)
    end
    function Model(V::Lar.Points)
        Lar.Model(V, true)
    end
    function Model()
        nothing
    end
end

#-------------------------------------------------------------------------------
#   BASIC PROPERTIES
#-------------------------------------------------------------------------------

length(m::Lar.Model)        = size(m.G, 1)
size(m::Lar.Model)		    = size(m.G)
size(m::Lar.Model, i::Int)  = size(m.G, i)
Base.copy(m::Lar.Model)     = Lar.Model(m.G, m.T)
Base.deepcopy(m::Lar.Model) = Lar.Model(Base.deepcopy(m.G), Base.deepcopy(m.T))
==(m1::Lar.Model, m2::Lar.Model) = m1.G == m2.G && m1.T == m2.T

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

    if !isempty(m.T[1])
        m.T[1] = [m.T[1] spzeros(Int8, m.T[1].m, 1)]
    end

    if !m.topological
        for i = 2 : length(m)
            if !isempty(m.T[i])
                m.T[i] = [m.T[i] spzeros(Int8, m.T[i].m, 1)]
            end
        end
    end

end

function deleteModelVertex!(m::Lar.Model, v::Int)::Nothing
    deleteModelVertices!(m, [v])
end

function deleteModelVertices!(m::Lar.Model, vs::Array{Int,1})::Nothing
    tokeepV = setdiff(collect(1 : size(m, 2)), vs);
    m.G     = m.G[:, tokeepV]

    if !isempty(m.T[1])
        # From EV must be deleted all the columns linked to a deleted vertex
        #  and all rows related to a dangling edge
        todelE = m.T[1][:, vs] != 0
        todelE = [l for l = 1 : todelE.m if sum(todelE[l, :]) != 0]
        m.T[1] = m.T[1][:, tokeepV]
        Lar.deleteModelEdges!(m, todelE)
    end

    if !m.topological
        for i = 2 : length(m)
            if !isempty(m.T[i])
                m.T[i] = m.T[i][:, tokeepV]
            end
        end
    end
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

function deleteModelEdges!(m::Lar.Model, es::Array{Int,1})::Nothing
    tokeepE = setdiff(collect(1 : m.T[1].m), es)
    m.T[1]  = m.T[1][es, :]
    #=
    EVs     = m.T[1]
    EVs     = EVs[  setdiff(collect(1 : EVs.m), ev),  :  ]
    tokeepV = [length(EVs[:,i].nzval) != 0  for i = 1 : EVs.n]
    m.G     = m.G[:, tokeepV]
    m.T[1]  = EVs[:, tokeepV]
    =#

    if m.topological && !isempty(m.T[2])
        m.T[2] = m.T[2][:, tokeepE]
    end

    return
end

#-------------------------------------------------------------------------------
#   MODEL MANIPULATION
#-------------------------------------------------------------------------------

function modelPurge!(m::Lar.Model, depth::Int = 0)::Nothing
    depth == 0 || throw(ArgumentError("NotCoded"))
    todel = [i for i = 1 : size(m, 2) if length(m.T[1][:, i].nzval) == 0]
    if !isempty(todel)
        Lar.deleteModelVertices!(m, todel)
    end
    return
end

#-------------------------------------------------------------------------------
#   MODEL CONVERSION TOPOLOGICAL ↔ GEOMETRICAL
#-------------------------------------------------------------------------------

function getGeometricalModel(tm::Lar.Model)::Lar.Model
    tm.topological || throw(ArgumentError("Model must be Topological"))
    gm = deepcopy(tm)
    gm.topological = false
    for i = 2 : length(gm)
        for j = i-1 : -1 : 1
            gm.T[i] = abs.(gm.T[i])*gm.T[j]
        end
    end
    return gm
end

function getTopologicalModel(gm::Lar.Model)::Lar.Model
    !gm.topological || throw(ArgumentError("Model must be Geometrical"))
    tm = deepcopy(gm)
    tm.topological = true
    for i = 2 : length(tm)
        for j = i-1 : -1 : 1
            #...
            throw(Error("Not Coded"))
        end
    end
    return tm
end
