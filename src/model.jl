import Base: +, length, size, ==, isempty
import ViewerGL
GL = ViewerGL

#-------------------------------------------------------------------------------
#   MODEL STRUCT DEFINITION
#-------------------------------------------------------------------------------

mutable struct Model
    G::Lar.Points
    T::Array{Lar.ChainOp, 1}

    function Model(V::Lar.Points, T::Array{Lar.ChainOp, 1})
        dim, npts = size(V)
        dim > 0 ||
            throw(ArgumentError("At least one point is needed."))
        length(T) == dim ||
            throw(ArgumentError("Topology is not coherent with Geometry."))

        isempty(T[1]) || size(T[1], 2) == npts ||
            throw(ArgumentError("Topology not coherent with Geometry."))
        for i = 2 : dim
            isempty(T[i-1]) || isempty(T[i]) ||
                size(T[i], 2) == size(T[i-1],1) ||
                throw(ArgumentError("Topology not coherent with Topology."))
        end

        new(V, T)
    end

    function Model(V::Lar.Points)
        T = convert(Array{Lar.ChainOp,1},
            [SparseArrays.spzeros(Int8, 0, 0) for i = 1 : size(V, 1)]
        )
        T[1] = convert(Lar.ChainOp, SparseArrays.spzeros(Int8, 0, size(V,2)))
        Lar.Model(V, T)
    end

    function Model()
        nothing
    end
end

#-------------------------------------------------------------------------------
#   BASIC PROPERTIES
#-------------------------------------------------------------------------------

length(m::Lar.Model)               = size(m.G, 1)
size(m::Lar.Model)		           = size(m.G)
size(m::Lar.Model, i::Int)         = i == 0 ? size(m.G   ) : size(m.T[i]   )
size(m::Lar.Model, i::Int, j::Int) = i == 0 ? size(m.G, j) : size(m.T[i], j)
==(m1::Lar.Model, m2::Lar.Model)   = m1.G == m2.G && m1.T == m2.T
isempty(m::Lar.Model, d::Int)      = isempty(m.T[d])

#-------------------------------------------------------------------------------
#   OTHER BASE REIMPLEMENTATIONS
#-------------------------------------------------------------------------------

Base.copy(m::Lar.Model)     = Lar.Model(m.G, m.T)
Base.deepcopy(m::Lar.Model) = Lar.Model(Base.deepcopy(m.G), Base.deepcopy(m.T))

#-------------------------------------------------------------------------------
#   GEOMETRY MANIPULATION
#-------------------------------------------------------------------------------

function addModelVertices!(m::Lar.Model, V::Lar.Points)::Nothing
    length(m) == size(V, 1) ||
        throw(ArgumentError("Point dimension mismatch."))
    m.G = [m.G V];

    m.T[1] = SparseArrays.sparse(
        findnz(m.T[1])...,
        size(m, 1, 1),
        size(m, 1, 2) + size(V, 2)
    )
    return
end

function addModelVertex!(m::Lar.Model, v::Array{Float64, 1})::Nothing
    return addModelVertices!(m, v[:, :])
end

function deleteModelVertex!(m::Lar.Model, v::Int)::Nothing
    deleteModelVertices!(m, [v])
end

function deleteModelVertices!(m::Lar.Model, vs::Array{Int,1})::Nothing
    tokeepV = setdiff(collect(1 : size(m, 0, 2)), vs);
    m.G     = m.G[:, tokeepV]

    if !isempty(m.T[1])
        # From EV must be deleted all the columns linked to a deleted vertex
        #  and all rows related to a dangling edge
        todelE = abs.(m.T[1][:, vs])
        todelE = [l for l = 1 : todelE.m if sum(todelE[l, :]) != 0]
        m.T[1] = m.T[1][:, tokeepV]
        if !isempty(todelE) Lar.deleteModelCells!(m, 1, todelE) end
    end

    return
end

#-------------------------------------------------------------------------------
#   TOPOLOGY MANIPULATION
#-------------------------------------------------------------------------------

function addModelCells!(m::Lar.Model, deg::Int, cs::Lar.ChainOp)::Nothing
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    (cs.n == size(m, deg - 1, 1) && deg > 1) ||
        (cs.n == size(m, deg - 1, 2) && deg == 1) ||
        throw(ArgumentError("Incoherent Chain Operator size"))

    # Adding input cells to the topological structure
    m.T[deg] = [m.T[deg]; cs]

    # Adding slot for the new cells to the higher order cells
    if deg < length(m)
        m.T[deg + 1] = [
            m.T[deg + 1] SparseArrays.spzeros(Int8, size(m, deg + 1, 1), cs.m)
        ]
    end

    return
end
function addModelCell!(m::Lar.Model, deg::Int, c::Lar.Cell)::Nothing
    Lar.addModelCells!(m, deg, convert(Lar.ChainOp, c))
end
function addModelCells!(m::Lar.Model, deg::Int, cs::Lar.Cells)::Nothing
    I = Array{Int,1}()
    J = Array{Int,1}()
    K = Array{Int8,1}()
    for i = 1 : length(cs)
        for j = 1 : length(cs[i])
            push!(I, i)
            push!(J, cs[i][j])
            push!(K, 1)
        end
    end

    scs = SparseArrays.sparse(I, J, K, length(cs), size(m, deg, 2));

    return addModelCells!(m, deg, scs);
end
function addModelCell!(m::Lar.Model, deg::Int, c::Array{Int64,1})::Nothing
    Lar.addModelCells!(m, deg, [c])
end


function deleteModelCells!(m::Lar.Model, deg::Int, cs::Array{Int, 1})::Nothing
    !isempty(cs) || return
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    !isempty(m.T[deg]) ||
        throw(ArgumentError("There are no cells of degree $deg"))
    max(cs...) ≤ size(m, deg, 1) ||
        throw(ArgumentError("There are not $(max(cs...)) cells of degree $deg"))

    # Removing `cs` rows from `m.T[deg]`
    tokeep     = setdiff(collect(1 : size(m, deg, 1)), cs)
    m.T[deg]   = m.T[deg][tokeep, :]

    if deg == length(m) return end

    # Removing `cs` cols from `m.T[deg+1]` by checking if some `deg+1` cell
    #  has to be removed too.
    # An higher ord cell must be removed if it contains a removed lower ord cell
    todelHo    = m.T[deg + 1][:, cs]
    m.T[deg+1] = m.T[deg + 1][:, tokeep]
    if isempty(m.T[deg+1]) return end
    todelHo    = [l for l = 1 : todel.m if sum(todelE[l, :]) != 0]
    if !isempty(todelHo) Lar.deleteModelCells!(m, deg+1, todelHo) end
    return
end

function deleteModelCell!(m::Lar.Model, deg::Int, c::Int)::Nothing
    Lar.deleteModelCells!(m, deg, [c])
end

function getModelLoCell(m::Lar.Model, deg::Int, c::Int)::Array{Int, 1}
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    return m.T[deg][c, :].nzind
end

function getModelCellVertices(m::Lar.Model, deg::Int, c::Int, ret_idx=false)
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    set = [c]
    for d = deg : -1 : 1
        set = ∪([m.T[d][el, :].nzind for el in set]...)
    end

    if ret_idx == true
        return (map(i -> m.G[:, set[i]], 1:length(set)), set)
    end
    return map(i -> m.G[:, set[i]], 1:length(set))
end

function getModelCellGeometry(m::Lar.Model, deg::Int, c::Int, ret_idx=false)
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    set = [c]
    for d = deg : -1 : 1
        set = ∪([m.T[d][el, :].nzind for el in set]...)
    end

    if ret_idx == true
        return (m.G[:, set], set)
    end
    return m.G[:, set]
end

#-------------------------------------------------------------------------------
#   MODEL MANIPULATION
#-------------------------------------------------------------------------------

"""
    modelPurge!(m::Lar.Model, [depth::Int = 0])::Nothing

Purges the model Geometry from unused elements.

Purges each Geometry element that is not related to a topological Cell.
If the optional argument `depth` is specified than it also purges all the
topological cells of degree lower or equal to depth that are not related
to an higher order topological cell.
"""
function modelPurge!(m::Lar.Model, depth::Int = 0)::Nothing
    depth < length(m) ||
        throw(ArgumentError("The model do not have degree $(depth+1)"))

    for d = depth : -1 : 1
        todel = [i for i = 1 : size(m, d+1, 2) if sum(abs.(m.T[d+1][:,i])) == 0]
        if !isempty(todel) Lar.deleteModelCells!(m, d, todel) end
    end
    todel = [i for i = 1 : size(m, 1, 2) if sum(abs.(m.T[1][:, i])) == 0]
    if !isempty(todel) Lar.deleteModelVertices!(m, todel) end
    return
end

#-------------------------------------------------------------------------------
#   MODEL CONVERSION TOPOLOGICAL ↔ GEOMETRICAL
#-------------------------------------------------------------------------------

#=
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
=#


#TODO check if removing sparse zeros is necessary



#-------------------------------------------------------------------------------
#   MODEL INTERACTION
#-------------------------------------------------------------------------------

function uniteModels!(model::Lar.Model, m2::Lar.Model)::Nothing

    length(model) == length(m2) ||
        throw(ArgumentError("ERROR: Inconsistent models dimension!"))

    Lar.addModelVertices!(model, m2.G)
    for d = 1 : length(model)
        shift = size(model, d, 2) - size(m2, d, 2)
        I, J, K = SparseArrays.findnz(m2.T[d])
        J = J .+ shift
        Lar.addModelCells!(model, d, SparseArrays.sparse(
            I, J, K, size(m2, d, 1), size(m2, d, 2) + shift
        ))
    end

    return
end

function uniteModels(m1::Lar.Model, m2::Lar.Model)::Lar.Model
    model = Lar.deepcopy(m1)
    Lar.uniteModels!(model, m2)
    return model
end

function mergeMultipleModels(models::Array{Lar.Model,1}; err=1e-6)::Lar.Model

    model = deepcopy(models[1]);

    for i = 2 : length(models)
        Lar.uniteModels!(model, models[i])
    end

    return Lar.mergeModelVertices(model, err=err)
end

function mergeModelVertices(model::Lar.Model; err=1e-4)
    V = convert(Lar.Points, model.G')
    EV = model.T[1]
    FE = model.T[2]
    vertsnum = size(V, 1)
    edgenum = size(EV, 1)
    facenum = size(FE, 1)
    newverts = zeros(Int, vertsnum)
    # KDTree constructor needs an explicit array of Float64
    V = Array{Float64,2}(V)
    W = convert(Lar.Points, LinearAlgebra.transpose(V))
    kdtree = KDTree(W)
	# remove vertices congruent to a single representative
    todelete = []
    i = 1
    for vi in 1:vertsnum
        if !(vi in todelete)
            nearvs = Lar.inrange(kdtree, V[vi, :], err)
            newverts[nearvs] .= i
            nearvs = setdiff(nearvs, vi)
            todelete = union(todelete, nearvs)
            i = i + 1
        end
    end
    nV = V[setdiff(collect(1:vertsnum), todelete), :]

    # translate edges to take congruence into account
    edges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
    oedges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
    for ei in 1:edgenum
        v1, v2 = EV[ei, :].nzind
        edges[ei] = Tuple{Int, Int}(sort([newverts[v1], newverts[v2]]))
        oedges[ei] = Tuple{Int, Int}(sort([v1, v2]))
    end
    nedges = union(edges)
    # remove edges of zero length
    nedges = filter(t->t[1]!=t[2], nedges)
    nedgenum = length(nedges)
    nEV = spzeros(Int8, nedgenum, size(nV, 1))

    etuple2idx = Dict{Tuple{Int, Int}, Int}()
    for ei in 1:nedgenum
    	begin
        	nEV[ei, collect(nedges[ei])] .= 1
        	nEV
        end
        etuple2idx[nedges[ei]] = ei
    end
    for e in 1:nedgenum
    	v1,v2 = findnz(nEV[e,:])[1]
    	nEV[e,v1] = -1; nEV[e,v2] = 1
    end

    # compute new faces to take congruence into account
    faces = [[
        map(x->newverts[x], FE[fi, ei] > 0 ? oedges[ei] : reverse(oedges[ei]))
        for ei in FE[fi, :].nzind
    ] for fi in 1:facenum]


    visited = []
    function filter_fn(face)

        verts = []
        map(e->verts = union(verts, collect(e)), face)
        verts = Set(verts)

        if !(verts in visited)
            push!(visited, verts)
            return true
        end
        return false
    end

    nfaces = filter(filter_fn, faces)

    nfacenum = length(nfaces)
    nFE = spzeros(Int8, nfacenum, size(nEV, 1))

    for fi in 1:nfacenum
        for edge in nfaces[fi]
            ei = etuple2idx[Tuple{Int, Int}(sort(collect(edge)))]
            nFE[fi, ei] = sign(edge[2] - edge[1])
        end
    end

    model = Lar.Model(convert(Lar.Points, nV'))
    Lar.addModelCells!(model, 1, nEV)
    Lar.addModelCells!(model, 2, nFE)
    return model
end

#=
function mergeModelVertices(model::Lar.Model; err=1e-6)
    dims = length(model)
    Gfloat = convert(Array{Float64,2}, model.G)
    kdtree = NearestNeighbors.KDTree(Gfloat)
    todel  = []
    tokeep = []
    I      = [Array{Int, 1}() for i = 1 : dims]
    J      = [Array{Int, 1}() for i = 1 : dims]
    K      = [Array{Int8,1}() for i = 1 : dims]

    # Collide Vertices and update 1-Cochain
    let nidx = 1;  for vidx = 1 : size(model, 0, 2)
        if !(vidx in todel)
            nearvs = NearestNeighbors.inrange(kdtree, Gfloat[:, vidx], err)
            # for d = 1 : dims
            if !isempty(model, 1)
                # @assert sum(sum(abs.(model.T[1][:, nearvs]), dims=2).>1) == 0 "IMPRECISION"
                newVL = sum(model.T[1][:, nearvs], dims=2)
                for i = 1 : length(newVL)  if abs(newVL[i]) == 1
                    push!(I[1], i)
                    push!(J[1], nidx)
                    push!(K[1], newVL[i])
                end  end
            end
            # end
            nidx = nidx + 1
            push!(todel, setdiff(nearvs, vidx)...)
            push!(tokeep, vidx)
        end
    end  end

    # Copy Upper Degree Topology
    for d = 2 : dims
        (I[d], J[d], K[d]) = SparseArrays.findnz(model.T[d])
    end

    model = Lar.Model(model.G[:, tokeep])


    # Collide Topology                                                          #TODO check
    let lolength = length(tokeep);  for d = 1 : dims  if !isempty(I[d])
        # Define compact indices for `d` dimension

        # Recostructiong Topology of degree `d`
        Td = SparseArrays.sparse(I[d], J[d], K[d], max(I[d]...), lolength)
        # Convert Topology to Array of Arrays
        aaT = [SparseArrays.findnz(Td[row, :])[1] for row = 1 : size(Td, 1)]
        nEl = length(aaT)
        # Find new Topology as array of arrays
        newT = sort(unique(aaT))
        # get the indices of the ordered topology
        ordT = sortperm([aaT; newT])
        # building the new Topology as SparseMatrix
        In = Array{Int , 1}()
        Jn = Array{Int , 1}()
        Kn = Array{Int8, 1}()
        let newidx = 1;  for eidx = 1 : length(ordT)  if newidx == 1 || ordT[eidx] < ordT[eidx-1]
            #                            v---co_edge------
            Je, Ke = SparseArrays.findnz(Td[ordT[eidx], :])
            push!(In, (ones(Int, length(Je)) * newidx)...)
            push!(Jn, Je...)
            push!(Kn, Ke...)
            newidx = newidx+1
        end  end  end
        Tn = SparseArrays.sparse(In, Jn, Kn, length(newT), lolength)
        # building the map originalEdge -> (collidedEdge, orientation)
        perm = [(0, 0) for i = 1 : nEl]
        let newidx = 1;  for e in ordT
            e <= nEl ? perm[e] = (newidx, Td[e, :] == Tn[newidx, :] ? 1 : -1) : newidx = newidx + 1
        end  end

        # updating HO topology indices
        for i = 1 : length(J[d+1])
            J[d+1][i], K[d+1][i] = perm[J[d+1][i]] .* (1, K[d+1][i])
        end

        # setting the new topology in the model
        Lar.addModelCells!(model, d, Tn);

        # getting the next dimension topology correct size
        lolength = size(model, d, 1)
    end  end  end

    # Last checks
    for d = 1 : length(model)
        todel = Array{Int, 1}()
        for i = 1 : size(model, d, 1)
            if sum(abs.(model.T[d][i, :])) <= d
                push!(todel, i)
            end
        end
        Lar.deleteModelCells!(model, d, todel)
    end

    return model
end
=#

#-------------------------------------------------------------------------------
#   MODEL REPRESENTATION
#-------------------------------------------------------------------------------

function view2DModel(model::Lar.Model, exp::Float64 = 1.)
    # visualization of numbered arrangement
    #VV = [[k] for k = 1 : size(model, 0, 2)]
    #EV = LAR.cop2lar(model.T[1]);
    #FE = LAR.cop2lar(model.T[2]);
    #FV = cat([[EV[e] for e in face] for face in FE])

    # final solid visualization
    FVs = Lar.triangulate2D(model)
    GL.VIEW(GL.GLExplode(model.G,FVs,exp,exp,exp,99,1));

    # polygonal face boundaries
    EVs = Lar.FV2EVs(model)
    GL.VIEW(GL.GLExplode(model.G,EVs,exp,exp,exp,3,1));
end

function viewModel(model::Lar.Model, depth::Int64 = -1, exp::Float64 = 1.)
    if depth == -1  depth = length(model)  end

    if depth == 0
        GL.VIEW( [GL.GLPoints(convert(Lar.Points, model.G'))] )
    end

    if depth == 1
        if length(model) < 3 || !isempty(model.T[3])
            copCE = abs.(model.T[3]) * abs.(model.T[2]) .!= 0
            comp = [(copCE[c, :]).nzind for c = 1 : size(copCE, 1)]
        elseif !isempty(model.T[2])
            comp = [(model.T[2][c, :]).nzind for c = 1 : size(model.T[2], 1)]
        else
            comp = Lar.Arrangement.biconnected_components(model.T[1])
        end

        mesh = [
            GL.GLGrid(
                model.G,
                LAR.cop2lar(model.T[1][comp[i], :]),
                GL.COLORS[(i-1)%12+1],
                1.
            )
            for i = 1 : length(comp)
        ]
        GL.VIEW( mesh )
    end

    if depth == 2
        model = deepcopy(model)
        modelPurge!(model, 1)
        FVs = Lar.triangulate2D(model)
        GL.VIEW(GL.GLExplode(model.G,FVs,exp,exp,exp,99,1));
    end

    if depth == 3
        GL.VIEW(  )
    end
end

#-------------------------------------------------------------------------------
#   MODEL EXPORT
#-------------------------------------------------------------------------------

function jlexportModel(
        model::Lar.Model,
        filename::String,
        modelname::String = "model"
    )::Nothing

    open(filename, "w") do f
        # Write imports
        write(f, "# using Lar\n")
        write(f, "using SparseArrays\n\n")

        # Write Geometry
        write(f, "#=== GEOMETRY ===#\n")
        write(f, "$modelname = Lar.Model([")
        for d = 1 : length(model)
            write(f, "\n\t")
            for p in model.G[d, :]  write(f, "$p ")  end
        end
        write(f, "\n]);\n\n")

        # Write Topology
        for d = 1 : length(model)
            write(f, "#=== TOPOLOGY : $d-CELLS ===#\n")
            I, J, K = findnz(model.T[d])
            write(f, "I = Array{Int64}([")
            for i in I  write(f, "$i, ")  end
            write(f, "]);\nJ = Array{Int64}([")
            for j in J  write(f, "$j, ")  end
            write(f, "]);\nK = Array{Int8}([")
            for k in K  write(f, "$k, ")  end
            write(f, "]);\n\n")
            write(f, "Lar.addModelCells!($modelname, $d, SparseArrays.sparse(")
            write(f, "I, J, K, $(size(model, d, 1)), $(size(model, d, 2))")
            write(f, "))\n\n\n")
        end
    end

    return
end
