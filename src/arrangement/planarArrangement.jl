#using Lar
#Lar = Lar;
LarA = Lar.Arrangement;

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 1 METHODS
#    frag_edge_channel
#    frag_edge
#    intersect_edges
#    merge_vertices!
#-------------------------------------------------------------------------------

"""                                                                             ## UNTESTED
    frag_edge_channel(in_chan, out_chan, model::Lar.Model, bigPI)

Utility function for parallel edge fragmentation.

This function handles the edge fragmentation in the first part of arrangement's
algorithmic pipeline (see also [`Lar.Arrangement.planar_arrangement_1`](@ref))
if multiprocessing computation has been enabled.
In order to do so, it needs two `Distributed.RemoteChannel`s, one with
the inputs and one for outputs.

See also: [`Lar.Arrangement.planar_arrangement_1`](@ref)
It uses: [`Lar.Arrangement.frag_edge`](@ref)
"""
function frag_edge_channel(
        in_chan::Distributed.RemoteChannel{Channel{Int64}},
        out_chan::Distributed.RemoteChannel{Channel{Tuple}},
        model::Lar.Model,
        bigPI::Array{Array{Int64,1},1}
    )
    run_loop = true
    while run_loop
        edgenum = take!(in_chan)
        if edgenum != -1
            put!(out_chan, (edgenum, frag_edge(model, edgenum, bigPI)))
        else
            run_loop = false
        end
    end
end

"""
    frag_edge(model::Lar.Model, edge_idx::Int, bigPI::Array)

Splits the `edge_idx`-th edge of the Edge Topology `model.T[1]`.

This method splits the `edge_idx`-th edge between the Edges in the Topology
`EV = model.T[1]` into several parts by confronting it with the others that
intersect its bounding box (see also [`Lar.Arrangement.intersect_edges`](@ref)).
Returns a set made of the detected vertices over the edge (with redundancies)
and the associated cochain (without redundancies).

See also:
 - [`Lar.Arrangement.planar_arrangement_1`](@ref)
 - [`Lar.Arrangement.frag_edge_channel`](@ref)
It uses: [`Lar.Arrangement.intersect_edges`](@ref)
---
# Examples
```jldoctest
julia> model = Lar.Model(hcat([
    [1.0, 0.0], [0.0, 1.0], [0.0, 0.5], [0.5, 1.0], [1.0, 1.0]
]...));
julia> model.T[1] = Lar.coboundary_0([[1, 2], [2, 5], [3, 4], [4, 5]]);
julia> bigPI = Lar.spaceindex((V, EV));
julia> LarA.frag_edge(model, 1, bigPI)[1]
4×2 Array{Float64,2}:
 1.0   0.0
 0.0   1.0
 0.25  0.75
 0.0   1.0
julia> LarA.frag_edge(model, 1, bigPI)[2]
2×4 SparseMatrixCSC{Int8,Int64} with 4 stored entries:
  [1, 1]  =  1
  [2, 2]  =  1
  [1, 3]  =  1
  [2, 3]  =  1
```
"""
function frag_edge(
        model::Lar.Model,
        edge_idx::Int,
        bigPI::Array{Array{Int64,1},1}
    )::Tuple{Lar.Points, Lar.ChainOp}
    V = convert(Lar.Points, model.G')
    EV = model.T[1]
    alphas = Dict{Float64, Int}()
    edge = EV[edge_idx, :]
    verts = V[edge.nzind, :]
    for i in bigPI[edge_idx]
        if i != edge_idx # && edge_idx in bigPI[i]
            intersection = LarA.intersect_edges(V, edge, EV[i, :])
            for (point, alpha) in intersection
                verts = [verts; point]
                alphas[alpha] = size(verts, 1)
            end
        end
    end
    alphas[0.0], alphas[1.0] = [1, 2]
    alphas_keys = sort(collect(keys(alphas)))
    edge_num = length(alphas_keys)-1
    verts_num = size(verts, 1)
    ev = SparseArrays.spzeros(Int8, edge_num, verts_num)
    for i in 1:edge_num
        ev[i, alphas[alphas_keys[i]]] = 1
        ev[i, alphas[alphas_keys[i+1]]] = 1
    end
    return verts, ev
end


"""
    intersect_edges(V::Lar.Points, edge1::Lar.Cell, edge2::Lar.Cell)

Finds the intersection points (if there exists) between the two given edges.

Compute the new points over `edge1` given by the intersection with `edge2`.
Note that vertices of `edge1` are therefore never detected (see example 2).

See also: [`Lar.Arrangement.frag_edge`](@ref)
---
# Examples
```jldoctest
# Cross
julia> V = [1.0 0.0; 0.0 1.0; 0.0 0.5; 0.5 1.0];
julia> copEV = Lar.coboundary_0([[1, 2], [3, 4]]);
julia> LarA.intersect_edges(V, copEV[1, :], copEV[2, :])
1-element Array{Tuple{Array{T,2} where T,Float64},1}:
 ([0.25 0.75], 0.75)
```

```jldoctest
# Collinear
julia> V = [1.0 0.0; 0.0 1.0; 0.75 0.25; 0.5 0.5];
julia> copEV = Lar.coboundary_0([[1, 2], [3, 4]]);
julia> Lar.Arrangement.intersect_edges(V, copEV[1, :], copEV[2, :])
2-element Array{Tuple{Array{T,2} where T,Float64},1}:
 ([0.75 0.25], 0.25)
 ([0.5 0.5], 0.5)
julia> Lar.Arrangement.intersect_edges(V, copEV[2, :], copEV[1, :])
0-element Array{Tuple{Array{T,2} where T,Float64},1}
```
"""
function intersect_edges(V::Lar.Points, edge1::Lar.Cell, edge2::Lar.Cell)
    err = 10e-8

    x1, y1, x2, y2 = vcat(map(c->V[c, :], edge1.nzind)...)
    x3, y3, x4, y4 = vcat(map(c->V[c, :], edge2.nzind)...)
    ret = Array{Tuple{Lar.Points, Float64}, 1}()

    v1 = [x2-x1, y2-y1];
    v2 = [x4-x3, y4-y3];
    v3 = [x3-x1, y3-y1];
    ang1 = dot(LA.normalize(v1), LA.normalize(v2))
    ang2 = dot(LA.normalize(v1), LA.normalize(v3))
    parallel = 1-err < abs(ang1) < 1+err
    colinear = parallel && (1-err < abs(ang2) < 1+err || -err < norm(v3) < err)
    if colinear
        o = [x1 y1]
        v = [x2 y2] - o
        alpha = 1/dot(v,v')
        ps = [x3 y3; x4 y4]
        for i in 1:2
            a = alpha*dot(v',(reshape(ps[i, :], 1, 2)-o))
            if 0 < a < 1
                push!(ret, (ps[i:i, :], a))
            end
        end
    elseif !parallel
        denom = (v2[2])*(v1[1]) - (v2[1])*(v1[2])
        a = ((v2[1])*(-v3[2]) - (v2[2])*(-v3[1])) / denom
        b = ((v1[1])*(-v3[2]) - (v1[2])*(-v3[1])) / denom

        if -err < a < 1+err && -err <= b <= 1+err
            p = [(x1 + a*(x2-x1))  (y1 + a*(y2-y1))]
            push!(ret, (p, a))
        end
    end
    return ret
end


"""
merge_vertices!(
    model::Lar.model,
    [ edge_map::Array{Array{Int64,1},1} = [[-1]] ],
    [ err::Float64 = 1e-4 ]
)::Nothing

Modify `model` by compacting its points closer than `err` in a single one.

This method check one at time each vertex ``v`` in `model.G` and identifies
each other vertex within `err` with ``v`` itself.
The edges cochain `model.T[1]` is coherently modified (multiple edges between
two vertices are not allowed).
If an `edge_map` is given in input (this could be usefull during the planar
arrangement), then also the map is coherently modified and given back in output.

See also: [`Lar.Arrangement.planar_arrangement_1`](@ref)
---
# Examples
```jldoctest
julia> model = Lar.Model([
    0.5 0.0 0.5 1.0 0.5 1.0
    0.5 0.0 0.5 1.0 0.5 1.0
]);
julia> model.T[1] = Lar.coboundary_0([[1, 4], [3, 2], [5, 6], [1, 6], [5, 3]]);
julia> LarA.merge_vertices!(model, [[-1]])

julia> model.G
2×3 Array{Float64,2}:
 0.5  0.0  1.0
 0.5  0.0  1.0

julia> model.T[1]
2×3 SparseMatrixCSC{Int8,Int64} with 4 stored entries:
  [1, 1]  =  1
  [2, 1]  =  1
  [2, 2]  =  1
  [1, 3]  =  1
```
"""
function merge_vertices!(
        model::Lar.Model,
        edge_map::Array{Array{Int64,1},1} = [[-1]],
        err::Float64 = 1e-4
    )::Nothing
    V = convert(Lar.Points, model.G')
    EV = model.T[1]
    vertsnum = size(V, 1)
    edgenum = size(EV, 1)
    newverts = zeros(Int, vertsnum)
    # KDTree constructor needs an explicit array of Float64
    V = Array{Float64,2}(V)
    kdtree = KDTree(permutedims(V))

    # merge congruent vertices
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

    # merge congruent edges
    edges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
    oedges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
    for ei in 1:edgenum
        v1, v2 = EV[ei, :].nzind
        edges[ei] = Tuple{Int, Int}(sort([newverts[v1], newverts[v2]]))
        oedges[ei] = Tuple{Int, Int}(sort([v1, v2]))
    end
    nedges = union(edges)
    nedges = filter(t->t[1]!=t[2], nedges)
    nedgenum = length(nedges)
    nEV = spzeros(Int8, nedgenum, size(nV, 1))
    # maps pairs of vertex indices to edge index
    etuple2idx = Dict{Tuple{Int, Int}, Int}()
    # builds `edge_map`
    for ei in 1:nedgenum
        nEV[ei, collect(nedges[ei])] .= 1
        etuple2idx[nedges[ei]] = ei
    end
    if edge_map != [[-1]]
        for i in 1:length(edge_map)
            row = edge_map[i]
            row = map(x->edges[x], row)
            row = filter(t->t[1]!=t[2], row)
            row = map(x->etuple2idx[x], row)
            edge_map[i] = row
        end
    end
    # return new vertices and new edges
    model.G = convert(Lar.Points, nV')
    model.T[1] = nEV
    return
end


#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - CLEAN DECOMPOSITION
#-------------------------------------------------------------------------------

"""
    cleandecomposition(model, sigma, edge_map)

                                                                                ## TODO

This function clears the `model` from all edges outside ``σ``.

This function takes a `Lar.Chain` ``σ`` and drops all the edges from
the topology `model.T[1]` that are outside from ``σ``.
gives back the dropped model and and a vector `todel` of the 2-cells to drop.

See also: [`Lar.planar_arrangement`](@ref)
---
# Examples
```jldoctest
# Nested Triangles
julia> V = [
		0.0 0.0; 2.0 0.0; 4.0 0.0;
		1.0 1.5; 3.0 1.5; 2.0 3.0;
		2.0 -3.; -2. 3.0; 6.0 3.0
	];

julia> copEV = SparseArrays.sparse(Array{Int8, 2}([
    [1 0 1 0 0 0 0 0 0] #1 -> 1,3
    [1 0 0 0 0 1 0 0 0] #2 -> 1,6
    [0 0 1 0 0 1 0 0 0] #3 -> 3,6
    [0 1 0 1 0 0 0 0 0] #4 -> 2,4
    [0 1 0 0 1 0 0 0 0] #5 -> 2,5
    [0 0 0 1 1 0 0 0 0] #6 -> 4,5
    [0 0 0 0 0 0 1 1 0] #7 -> 7,8
    [0 0 0 0 0 0 1 0 1] #8 -> 7,9
    [0 0 0 0 0 0 0 1 1] #9 -> 8,9
    ]));

julia> σ = SparseArrays.sparse([0; 0; 0; 1; 1; 1; 0; 0; 0]);

todel, V, copEV = LarA.cleandecomposition(V, copEV, convert(Lar.Chain, σ))

Plasm.view(convert(Lar.Points, V'), Lar.cop2lar(copEV));
```
"""
function cleandecomposition(V, copEV, sigma, edge_map)
    # Deletes edges outside sigma area
    todel = []
    new_edges = []
    map(i->new_edges=union(new_edges, edge_map[i]), sigma.nzind)
    ev = copEV[new_edges, :]
    for e in 1:copEV.m
        if !(e in new_edges)

            vidxs = copEV[e, :].nzind
            v1, v2 = map(i->V[vidxs[i], :], [1,2])
            centroid = .5*(v1 + v2)

            if ! Lar.point_in_face(centroid, V, ev)
                push!(todel, e)
            end
        end
    end

    for i in reverse(todel)
        for row in edge_map

            filter!(x->x!=i, row)

            for j in 1:length(row)
                if row[j] > i
                    row[j] -= 1
                end
            end
        end
    end

    V, copEV = Lar.delete_edges(todel, V, copEV)
	return V,copEV
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - BICONNECTED COMPONENTS
#-------------------------------------------------------------------------------

"""
    biconnected_components(EV::Lar.ChainOp)

Compute sets of 2-cells on the same 3-cells biconnected components.

The method evaluates the 2-cells ``σ_i`` wich lies on the same 3-cells
biconnected component ``χ_j`` and gives back an array that contains an array
of ``σ_i`` for each ``j``.
Do note that if the cochain `EV` contains more copies of the same 2-cell then
it will be considered like a 3-cell.

Intermediate part of Planar Arrangement's algorithmic pipeline.

See also:
 - [`Lar.planar_arrangement`](@ref) for the complete pipeline.
 - [`Lar.Arragnement.planar_arrangement_2`](@ref).
---

# Examples
```jldoctest
julia> copEV = SparseArrays.sparse(Array{Int8, 2}([
        [1 1 0 0 0 0] #1 -> 1,2  |
        [1 0 1 0 0 0] #2 -> 1,3  |
        [1 0 0 1 0 0] #3 -> 1,4   |
        [1 0 0 0 1 0] #4 -> 1,5   |
        [1 0 0 0 0 1] #5 -> 1,6
        [0 1 1 0 0 0] #6 -> 2,3  |
        [0 0 0 1 1 0] #7 -> 4,5   |
    ]));
julia> LarA.biconnected_components(copEV)
2-element Array{Array{Int64,1},1}:
 [2, 6, 1]
 [4, 7, 3]
```

```jldoctest
julia> copEV = SparseArrays.sparse(Array{Int8, 2}([
        [1 1 0] #1 -> 1,2  |
        [1 1 0] #2 -> 1,2  |
        [1 0 1] #3 -> 1,2
    ]));
julia> LarA.biconnected_components(copEV)
1-element Array{Array{Int64,1},1}:
 [2, 1]
```
"""
function biconnected_components(EV::Lar.ChainOp)

    ps = Array{Tuple{Int, Int, Int}, 1}()
    es = Array{Tuple{Int, Int}, 1}()
    todel = Array{Int, 1}()
    visited = Array{Int, 1}()
    bicon_comps = Array{Array{Int, 1}, 1}()
    hivtx = 1

    function an_edge(point) # TODO: fix bug
        # error? : BoundsError: attempt to access 0×0 SparseMatrix ...
        edges = setdiff(EV[:, point].nzind, todel)
        if length(edges) == 0
            edges = [false]
        end
        edges[1]
    end

    function get_head(edge, tail)
        setdiff(EV[edge, :].nzind, [tail])[1]
    end

    function v_to_vi(v)
        i = findfirst(t->t[1]==v, ps)
        # seems findfirst changed from 0 to Nothing
        if typeof(i) == Nothing
            return false
        elseif i == 0
            return false
        else
            return ps[i][2]
        end
    end

    push!(ps, (1,1,1))
    push!(visited, 1)
    exit = false
    while !exit
        edge = an_edge(ps[end][1])
        if edge != false
            tail = ps[end][2]
            head = get_head(edge, ps[end][1])
            hi = v_to_vi(head)
            if hi == false
                hivtx += 1
                push!(ps, (head, hivtx, ps[end][2]))
                push!(visited, head)
            else
                if hi < ps[end][3]
                    ps[end] = (ps[end][1], ps[end][2], hi)
                end
            end
            push!(es, (edge, tail))
            push!(todel, edge)
        else
            if length(ps) == 1
                found = false
                pop!(ps)
                for i in 1:size(EV,2)
                    if !(i in visited)
                        hivtx = 1
                        push!(ps, (i, hivtx, 1))
                        push!(visited, i)
                        found = true
                        break
                    end
                end
                if !found
                    exit = true
                end

            else
                if ps[end][3] == ps[end-1][2]
                    edges = Array{Int, 1}()
                    while true
                        edge, tail = pop!(es)
                        push!(edges, edge)
                        if tail == ps[end][3]
                            if length(edges) > 1
                                push!(bicon_comps, edges)
                            end
                            break
                        end
                    end

                else
                    if ps[end-1][3] > ps[end][3]
                        ps[end-1] = (ps[end-1][1], ps[end-1][2], ps[end][3])
                    end
                end
                pop!(ps)
            end
        end
    end
    bicon_comps = sort(bicon_comps, lt=(x,y)->length(x)>length(y))
    return bicon_comps
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2 - COMPONENT GRAPH (TGW) METHODS
#    get_external_cycle
#    pre_containment_test
#    prune_containment_graph
#    transitive_reduction!
#-------------------------------------------------------------------------------

"""
    get_external_cycle(model::Lar.Model)

Evaluates the index of the external 3D-cell.

This method looks for and retrieve the external cycle of the `model`.
If the cell does not exist then it return a `Nothing` element.

See also: [`Lar.Arrangement.componentgraph`](@ref).
---
# Examples
```jldoctest
# ``K^4`` graph
julia> model = Lar.Model(hcat([[0.0,0.0], [4.0,0.0], [2.0,3.0], [2.0,1.5]]...));
julia> model.T[1] = SparseArrays.sparse(Array{Int8, 2}([
		[1 1 0 0] #1 -> 1,2
		[0 1 1 0] #2 -> 2,3
		[1 0 1 0] #3 -> 3,1
		[1 0 0 1] #4 -> 1,2
		[0 1 0 1] #5 -> 2,3
		[0 0 1 1] #6 -> 3,1
	]));
julia> model.T[2] = SparseArrays.sparse(Array{Int8, 2}([
		[1 0 0 1 1 0] #1 -> 1,4,5
		[0 1 0 0 1 1] #2 -> 2,5,6
		[0 0 1 1 0 1] #3 -> 3,4,6
		[1 1 1 0 0 0] #4 -> 1,2,3 External
	]));

julia> LarA.get_external_cycle(model)
4

julia> model.T[2] = model.T[2][1:3,:];
julia> typeof(LarA.get_external_cycle(model))
Nothing
```
"""
function get_external_cycle(V::Lar.Points, EV::Lar.ChainOp, FE::Lar.ChainOp)::Int64
    FV = abs.(FE)*EV
    vs = sparsevec(mapslices(sum, abs.(EV), dims=1)').nzind
    minv_x1 = maxv_x1 = minv_x2 = maxv_x2 = pop!(vs)
    for i in vs
        if V[i, 1] > V[maxv_x1, 1]
            maxv_x1 = i
        elseif V[i, 1] < V[minv_x1, 1]
            minv_x1 = i
        end
        if V[i, 2] > V[maxv_x2, 2]
            maxv_x2 = i
        elseif V[i, 2] < V[minv_x2, 2]
            minv_x2 = i
        end
    end
    cells = intersect(
        FV[:, minv_x1].nzind,
        FV[:, maxv_x1].nzind,
        FV[:, minv_x2].nzind,
        FV[:, maxv_x2].nzind
    )
    if length(cells) == 1
        return cells[1]
    else
        for c in cells
            if Lar.face_area(V, EV, FE[c, :]) < 0
                return c
            end
        end
    end
end

"""
    pre_containment_test(
        bboxes::Array{Tuple{Array{Float64,2},Array{Float64,2}},1}
    )

Generate the containment graph associated to `bboxes`.

The function evaluate a `SparseArrays{Int8, n, n}` PRE-containmnet graph meaning
it stores `1` when the bbox associated to the `i`-th component contains entirely
the bbox associated to the `j`-th one. (where `n = length(bboxes)`)

See also: [`Lar.Arrangement.componentgraph`](@ref).

---

# Examples
```jldoctest
# Planar Tiles
julia> bboxes = [
		([0.0 0.0], [1.0 1.0])
		([0.0 0.0], [0.4 0.4])
		([0.6 0.0], [1.0 0.4])
		([0.0 0.6], [0.4 1.0])
		([0.6 0.6], [1.0 1.0])
		([2.0 3.0], [2.0 3.0]) # no intersection
	];

julia> LarA.pre_containment_test(bboxes)
6x6 SparseMatrixCSC{Int8,Int64} with 4 stored entries:
  [2, 1]  =  1
  [3, 1]  =  1
  [4, 1]  =  1
  [5, 1]  =  1
```
"""
function pre_containment_test(
        bboxes#::Array{Tuple{Array{Float64,2},Array{Float64,2}},1}              #
    )
    n = length(bboxes)
    containment_graph = spzeros(Int8, n, n)

    for i in 1:n
        for j in 1:n
            if i != j && Lar.bbox_contains(bboxes[j], bboxes[i])
                containment_graph[i, j] = 1
            end
        end
    end

    return containment_graph
end

"""
    prune_containment_graph(n, V, EVs, shells, graph)

Prunes the containment `graph` from the non-included faces.

This method prunes the containment graph eliminating the included bouning boxes
that not corresponds to included faces.

Do note that this method expects to work on biconnectet clusters of faces.

See also: [`Lar.Arrangement.componentgraph`](@ref).
```
"""
function prune_containment_graph(n, V, EVs, shells, graph)

    for i in 1:n
        an_edge = shells[i].nzind[1]
        origin_index = EVs[i][an_edge, :].nzind[1]
        origin = V[origin_index, :]

        for j in 1:n
            if i != j
                if graph[i, j] == 1
                    shell_edge_indexes = shells[j].nzind
                    ev = EVs[j][shell_edge_indexes, :]

                    if !Lar.point_in_face(origin, V, ev)
                        graph[i, j] = 0
                    end
                end
             end
         end

     end
     return graph
end

"""
    transitive_reduction!(graph)

Prune the redundancies in containment `graph`.

The ``m``-layer containment `graph` is reduced to the corresponding 1-layer
containment `graph` where only the first level nested components are considered.

See also: [`Lar.Arrangement.componentgraph`](@ref).

---

# Examples
```jldoctest
julia> containment = [
		0 0 0 0 0 0
		0 0 0 0 0 0
		1 0 0 0 0 0
		1 0 0 0 0 0
		1 0 1 0 0 0
		1 0 1 0 1 0
	];

julia> LarA.transitive_reduction!(containment)

julia> containment
6×6 Array{Int64,2}:
 0  0  0  0  0  0
 0  0  0  0  0  0
 1  0  0  0  0  0
 1  0  0  0  0  0
 0  0  1  0  0  0
 0  0  0  0  1  0
```
"""
function transitive_reduction!(graph)
    n = size(graph, 1)
    for j in 1:n
        for i in 1:n
            if graph[i, j] > 0
                for k in 1:n
                    if graph[j, k] > 0
                        graph[i, k] = 0
                    end
                end
            end
        end
    end
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2 METHODS
#    componentgraph (aka TGW)
#    cell_merging
#-------------------------------------------------------------------------------

"""
    componentgraph(V, copEV, bicon_comps)

Topological Gift Wrapping algorithm on 2D skeletons.

This is the offline part of the TGW algorithm. It takes in input a model and its
biconnected components mapping and evaluates usefull informations:
 1. Number of biconnected components.
 2. Component Graph of the biconnected structure.
 3. The 1-cells structure (UNMODIFIED).                                         # Could be removed?
 4. Association between non-dangling 2-cells and their orientation (∀ component)
 5. Association between 3-cells and 2-cells (with orientation, ∀ component).
 6. Association between 3-cells and their orientation (∀ component).
 7. Shell bounding boxes of the components.

See also: [`Lar.Arrangement.planar_arrangement_2`](@ref).
It uses:
 - [`Lar.Arrangement.get_external_cycle`](@ref)
 - [`Lar.Arrangement.pre_containment_test`](@ref)
 - [`Lar.Arrangement.prune_containment_graph`](@ref)
 - [`Lar.Arrangement.transitive_reduction!`](@ref)
---
# Examples
```jldoctest
# Papillon
julia> V = [0.0 0.0; 0.0 3.0; 2.0 1.5; 4.0 0.0; 4.0 3.0];

julia> copEV = SparseArrays.sparse(Array{Int8, 2}([
		[1 1 0 0 0] #1 -> 1,2
		[0 1 1 0 0] #2 -> 2,3
		[1 0 1 0 0] #3 -> 3,1
		[0 0 1 1 0] #4 -> 3,4
		[0 0 0 1 1] #5 -> 4,5
		[0 0 1 0 1] #6 -> 3,5
	]));

julia> copFE = SparseArrays.sparse(Array{Int8, 2}([
		[1 1 1 0 0 0] #1 -> 1,2,3
		[0 0 0 1 1 1] #2 -> 4,5,6
		[1 1 1 1 1 1] #3 -> 1,2,3,4,5,6 External
	]));

julia> bicon_comps = [[1, 2, 3], [4, 5, 6]];

julia> LarA.componentgraph(V, copEV, bicon_comps)[1]
2

julia> LarA.componentgraph(V, copEV, bicon_comps)[2]
2×2 SparseMatrixCSC{Int8,Int64} with 0 stored entries

julia> LarA.componentgraph(V, copEV, bicon_comps)[4]
2-element Array{SparseMatrixCSC{Int8,Int64},1}:

  [1, 1]  =  -1
  [3, 1]  =  -1
  [1, 2]  =  1
  [2, 2]  =  -1
  [2, 3]  =  1
  [3, 3]  =  1

  [1, 3]  =  -1
  [3, 3]  =  -1
  [1, 4]  =  1
  [2, 4]  =  -1
  [2, 5]  =  1
  [3, 5]  =  1

julia> LarA.componentgraph(V, copEV, bicon_comps)[5]
2-element Array{SparseMatrixCSC{Int8,Int64},1}:

  [1, 1]  =  -1
  [1, 2]  =  -1
  [1, 3]  =  1

  [1, 1]  =  1
  [1, 2]  =  1
  [1, 3]  =  -1

julia> LarA.componentgraph(V, copEV, bicon_comps)[6]
2-element Array{SparseVector{Int8,Int64},1}:
   [1]  =  1
  [2]  =  1
  [3]  =  -1
   [1]  =  -1
  [2]  =  -1
  [3]  =  1

julia> LarA.componentgraph(V, copEV, bicon_comps)[7]
2-element Array{Any,1}:
 ([0.0 0.0], [2.0 3.0])
 ([2.0 0.0], [4.0 3.0])
```
"""
function componentgraph(V, copEV, bicon_comps)
    # arrangement of isolated components
	n = size(bicon_comps, 1)
   	shells = Array{Lar.Chain, 1}(undef, n)
	boundaries = Array{Lar.ChainOp, 1}(undef, n)
	EVs = Array{Lar.ChainOp, 1}(undef, n)
    # for each component
	for p in 1 : n
		ev = copEV[sort(bicon_comps[p]), :]
        # computation of 2-cells
		fe = Lar.Arrangement.minimal_2cycles(V, ev)
        # exterior cycle
		shell_num = Lar.Arrangement.get_external_cycle(V, ev, fe)
        # decompose each fe (co-boundary local to component)
		EVs[p] = ev
		tokeep = setdiff(1:fe.m, shell_num)
		boundaries[p] = fe[tokeep, :]
		shells[p] = fe[shell_num, :]
    end
    # computation of bounding boxes of isolated components
	shell_bboxes = []
	for i in 1:n
    	vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
   		push!(shell_bboxes, Lar.bbox(V[vs_indexes, :]))
	end
    # computation and reduction of containment graph
	containment_graph = pre_containment_test(shell_bboxes)
	containment_graph = prune_containment_graph(
        n, V, EVs, shells, containment_graph
    )
	transitive_reduction!(containment_graph)
	return n, containment_graph, V, EVs, boundaries, shells, shell_bboxes
end

"""
    cell_merging(n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)

Cells composing for the Topological Gift Wrapping algorithm.

This is the online part of the TGW algorithm.

See also: [`Lar.Arrangement.planar_arrangement_2`](@ref).
"""
function cell_merging(
        n, containment_graph, V, EVs, boundaries, shells, shell_bboxes
    )
    function bboxes(V::Lar.Points, indexes::Lar.ChainOp)
        boxes = Array{Tuple{Any, Any}}(undef, indexes.n)
        for i in 1:indexes.n
            v_inds = indexes[:, i].nzind
            boxes[i] = Lar.bbox(V[v_inds, :])
        end
        boxes
    end
    # initiolization
    sums = Array{Tuple{Int, Int, Int}}(undef, 0);
    # assembling child components with father components
    for father in 1:n
        if sum(containment_graph[:, father]) > 0
            father_bboxes = bboxes(
                V, abs.(EVs[father]')*abs.(boundaries[father]')
            )
            for child in 1:n
                if containment_graph[child, father] > 0
                    child_bbox = shell_bboxes[child]
                    for b in 1:length(father_bboxes)
                        if Lar.bbox_contains(father_bboxes[b], child_bbox)
                            push!(sums, (father, b, child))
                            break
                        end
                    end
                end
            end
        end
    end
    # offset assembly initialization
    EV = vcat(EVs...)
    edgenum = size(EV, 1)
    facenum = sum(map(x->size(x,1), boundaries))
    FE = spzeros(Int8, facenum, edgenum)
    shells2 = spzeros(Int8, length(shells), edgenum)
    r_offsets = [1]
    c_offset = 1
    # submatrices construction
    for i in 1:n
        min_row = r_offsets[end]
        max_row = r_offsets[end] + size(boundaries[i], 1) - 1
        min_col = c_offset
        max_col = c_offset + size(boundaries[i], 2) - 1
        FE[min_row:max_row, min_col:max_col] = boundaries[i]
        shells2[i, min_col:max_col] = shells[i]
        push!(r_offsets, max_row + 1)
        c_offset = max_col + 1
    end
    # offsetting assembly of component submatrices
    for (f, r, c) in sums
        FE[r_offsets[f]+r-1, :] += shells2[c, :]
    end

    return EV, FE
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 1
#-------------------------------------------------------------------------------


"""
	planar_arrangement_1(
        model::Lar.Model,
        [sigma],
        [return_edge_map],
        [multiproc]
    )

First part of arrangement's algorithmic pipeline.

This function computes the pairwise intersection between each edge of a given
2D cellular complex 1-skeleton. The computation is speeded up via the evaluation
of the Spatial Index. See [`Lar.spaceindex`](@ref).

See also: [`Lar.planar_arrangement`](@ref) for the complete pipeline.
It uses:
 - [`Lar.Arrangement.frag_edges_channel`](@ref)
 - [`Lar.Arrangement.frag_edges`](@ref)
 - [`Lar.Arrangement.merge_vertices!`](@ref)
---
# Examples
```jldoctest
julia> model = Lar.Model([
    0.0 0.5 0.0 0.5 0.3 1.0 0.3 1.0;
    0.0 0.0 1.0 1.0 0.5 0.5 1.0 1.0
]);
julia> model.T[1] = Lar.coboundary_0(
    [[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [5, 7], [6, 8]]
);
julia> (model, sigma, edge_map) =
    LarA.planar_arrangement_1(model, spzeros(Int8, 0), true);

julia> model.G
2×9 Array{Float64,2}:
 0.0  0.5  0.0  0.5  0.3  0.5  0.3  1.0  1.0
 0.0  0.0  1.0  1.0  1.0  0.5  0.5  0.5  1.0

julia> model.T[1]
11×9 SparseMatrixCSC{Int8,Int64} with 22 stored entries:
 [1 , 1]  =  1
 [4 , 1]  =  1
 ⋮
 [9 , 9]  =  1
 [11, 9]  =  1

julia> sigma
0-element SparseVector{Int8,Int64} with 0 stored entries

julia> edge_map
8-element Array{Array{Int64,1},1}:
 [1]
 [2, 3]
 ⋮
 [10]
 [11]

```
"""
function planar_arrangement_1(
        model::Lar.Model,
		sigma::Lar.Chain = spzeros(Int8, 0),
		return_edge_map::Bool = false,
		multiproc::Bool = false
    )
	# data structures initialization
    # V = convert(Lar.Points, model.G')
    # copEV = model.T[1]
	edgenum = size(model.T[1], 1)
	edge_map = Array{Array{Int, 1}, 1}(undef,edgenum)
	rV = Lar.Points(zeros(0, 2))
	rEV = SparseArrays.spzeros(Int8, 0, 0)
	finalcells_num = 0

	# spaceindex computation
	bigPI = Lar.spaceindex((model.G, Lar.cop2lar(model.T[1]))::Lar.LAR)

    # multiprocessing of edge fragmentation
    if (multiproc == true)
        in_chan = Distributed.RemoteChannel(()->Channel{Int64}(0))
        out_chan = Distributed.RemoteChannel(()->Channel{Tuple}(0))
        ordered_dict = DataStructures.SortedDict{Int64,Tuple}()
        @async begin
            for i in 1:edgenum
                put!(in_chan,i)
            end
            for p in Distributed.workers()
                put!(in_chan,-1)
            end
        end
        for p in Distributed.workers()
            @async Base.remote_do(
                frag_edge_channel, p, in_chan, out_chan, V, model.T[1], bigPI
            )
        end
        for i in 1:edgenum
            frag_done_job = take!(out_chan)
            ordered_dict[frag_done_job[1]] = frag_done_job[2]
        end
        for (dkey, dval) in ordered_dict
            i = dkey
            v, ev = dval
            newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
            edge_map[i] = newedges_nums
            finalcells_num += size(ev, 1)
            rV, rEV = Lar.skel_merge(rV, rEV, v, ev)
        end
    else
        # sequential (iterative) processing of edge fragmentation
        for i in 1:edgenum
            v, ev = LarA.frag_edge(model, i, bigPI)
            newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
            edge_map[i] = newedges_nums
            finalcells_num += size(ev, 1)
            rV = convert(Lar.Points, rV)
            rV, rEV = Lar.skel_merge(rV, rEV, v, ev)
        end
    end
    # merging of close vertices and edges (2D congruence)
    model = Lar.Model(convert(Lar.Points, rV'))
    model.T[1] = rEV
    LarA.merge_vertices!(model, edge_map)
    if !isnothing(sigma)
        sigmaN = sigma                                                          ## TODO
    end
    if return_edge_map
        return model, sigmaN, edge_map
    else
        return model, sigmaN
    end
end

"""
	planar_arrangement_2(model, [bicon_comps], [edge_map])

Second part of arrangement's algorithmic pipeline.

This function is the complete Topological Gift Wrapping (TGW) algorithm that is firstly
locally used in order to decompose the 2-cells and then globally to generate the 3-cells
of the arrangement of the ambient space ``E^3``.

During this process each dangling 2-cell is removed.
Do note that the isolated 1-cells are not removed by this procedure.

See also: [`Lar.planar_arrangement`](@ref) for the complete pipeline.
It uses:
 - [`Lar.Arrangement.biconnected_components`](@ref)
 - [`Lar.Arrangement.componentgraph`](@ref)
 - [`Lar.Arrangement.cell_merging`](@ref)
---

# Examples
```jldoctest
# Triforce
julia> model = Lar.Model(hcat([
    [0.0, 0.0],[2.0, 0.0],[4.0, 0.0],[1.0, 1.5],[3.0, 1.5],[2.0, 3.0],[3.0, 3.0]
    ]...));
julia> model.T[1] = SparseArrays.sparse(Array{Int8, 2}([
    [1 1 0 0 0 0 0] #1 -> 1,2
    [0 1 1 0 0 0 0] #2 -> 2,3
    [1 0 0 1 0 0 0] #3 -> 1,4
    [0 0 0 1 0 1 0] #4 -> 4,6
    [0 0 1 0 1 0 0] #5 -> 3,5
    [0 0 0 0 1 1 0] #6 -> 5,6
    [0 1 0 1 0 0 0] #7 -> 2,4
    [0 1 0 0 1 0 0] #8 -> 2,5
    [0 0 0 1 1 0 0] #9 -> 4,5
    [0 0 0 0 0 1 1]
    ]));

julia> bicon_comps = Lar.Arrangement.biconnected_components(model.T[1])
1-element Array{Array{Int64,1},1}:
 [8, 9, 7, 3, 4, 6, 5, 2, 1]

julia> V, EV, FE = Lar.Arrangement.planar_arrangement_2(model, bicon_comps);    ##Add numbering view
```
"""
function planar_arrangement_2(
        model::Lar.Model,
        bicon_comps,                                                            ##
        edge_map,                                                               ##
    )::Lar.Model

    edges = sort(union(bicon_comps...))
    todel = sort(setdiff(collect(1:size(model.T[1], 1)), edges))

    for i in reverse(todel)
        for row in edge_map

            filter!(x->x!=i, row)

            for j in 1:length(row)
                if row[j] > i
                    row[j] -= 1
                end
            end
        end
    end

    # test to remove
    @assert isequal(
        bicon_comps,
        Lar.Arrangement.biconnected_components(model.T[1])
    )

	# component graph
	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes =
        Lar.Arrangement.componentgraph(
            convert(Lar.Points, model.G'), model.T[1], bicon_comps
        )

	copEV, FE = Lar.Arrangement.cell_merging(
	   	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)

    model = Lar.Model(convert(Lar.Points, V'));
    model.T[1] = copEV;
    model.T[2] = FE;
	return model
end


#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE
#-------------------------------------------------------------------------------

"""
    planar_arrangement(
        model::Lar.Model,
        [sigma::Lar.Chain],
        [return_edge_map::Bool],
        [multiproc::Bool]
    )::Union{Lar.Model, Tuple{Lar.model, Array{Array{Int, 1}, 1}}}

Compute the arrangement on the given cellular complex 1-skeleton in 2D.

A cellular complex is arranged when the intersection of every possible pair
of cell of the complex is empty and the union of all the cells is the
whole Euclidean space.
The methods default input/output is a `Lar.Model` with at least Geometry and
1-cell Topoly.
Mutiprocessing is by default disabled.

It uses:
 - [`Lar.Arrangement.planar_arrangement_1`](@ref)
 - [`Lar.Arrangement.cleandecomposition`](@ref)
 - [`Lar.Arrangement.biconnected_components`](@ref)
 - [`Lar.Arrangement.planar_arrangement_2`](@ref)

## Additional arguments:
 - `sigma::Chain`: if specified, the method return only cells contained inside
        this cell. _Default_: empty cell.
 - `return_edge_map::Bool`: builds an `edge_map` which maps the input edges
        to the output's. _Defaults_: `false`.
 - `multiproc::Bool`: Runs the method in parallel mode. _Defaults_: `false`.
"""
function planar_arrangement(
        model::Lar.Model,
        sigma::Lar.Chain = spzeros(Int8, 0),                                    ## TEST
        return_edge_map::Bool = false,
        multiproc::Bool = false
    )::Union{Lar.Model, Tuple{Lar.Model, Array{Array{Int, 1}, 1}}}

    # Check dimensional condition
    @assert length(model) == 2;

    # Chek multiprocessing
    if multiproc && Distributed.nprocs() <= 1
        multiproc = false
        println("Setting Multiproc to False due to insufficient processors.")
    end

    #planar_arrangement_1
	model, sigma, edge_map =
        LarA.planar_arrangement_1(model, sigma, true, multiproc)

    # cleandecomposition
	if !isempty(sigma)
		model = Lar.Arrangement.cleandecomposition(model, sigma, edge_map)      ## TODO
	end

    bicon_comps = Lar.Arrangement.biconnected_components(model.T[1])
    # EV = Lar.cop2lar(copEV)
    # V,bicon_comps = Lar.biconnectedComponent((V,EV))

	if isempty(bicon_comps)
    	println("No biconnected components found.")
    	if (return_edge_map)
    	    return (model(), nothing);
    	else
    	    return model();
    	end
	end

    #Planar_arrangement_2
	model = Lar.Arrangement.planar_arrangement_2(
        model, bicon_comps, edge_map
    );

	if return_edge_map
	     return model, edge_map
	else
	     return model
	end
end
