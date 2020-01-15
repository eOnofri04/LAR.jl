using SparseArrays, LightXML, ViewerGL
GL = ViewerGL
include("./../src/Lar.jl")


#=== EXTRACT NODEs AND STREETs FROM OSM FILE ===#
xdoc  = parse_file("examples/map.osm")
xroot = root(xdoc)

node_list = get_elements_by_tagname(xroot, "node")
way_list  = get_elements_by_tagname(xroot, "way")

#=== EXTRACT VERTEX SET AND VERTEX POSITION MAP ===#

Vidxs = [parse(Int, attribute(v, "id")) for v in node_list]
Vmax  = max(Vidxs...)
Vmap  = SparseArrays.spzeros(Int, Vmax)
Vnum  = length(Vidxs)
for i = 1 : Vnum   Vmap[Vidxs[i]] = i   end


lats  = [parse(Float64, attribute(v, "lat")) for v in node_list]
lmin  = min(lats...)
lmax  = max(lats...)
lats  = (lats .- lmin) ./ (lmax - lmin)

lons  = [parse(Float64, attribute(v, "lon")) for v in node_list]
lmin  = min(lons...)
lmax  = max(lons...)
lons  = (lons .- lmin) ./ (lmax - lmin)

#=== BUILD LAR MODEL ===#

model = Lar.Model([lats'; lons'])

for street in way_list
    nodes = [
             Vmap[parse(Int, attribute(n, "ref"))]
             for n in get_elements_by_tagname(street, "nd")
            ]
    for i = 1 : length(nodes) - 1
        newcell = SparseArrays.spzeros(Int8, Vnum)
        newcell[nodes[  i]] = 1
        newcell[nodes[i+1]] = 1
        Lar.addModelCells!(model, 1, convert(Lar.ChainOp, newcell'))
    end
end


#===
    extraction of a box of
     0.31<lats<0.55
     0.65<lons<0.75
===#

toDelV = convert(Array{Int64,1},
    [i  for i = 1 : size(model, 0, 2)
        if model.G[1, i] < 0.31
        || model.G[1, i] > 0.35
        || model.G[2, i] < 0.65
        || model.G[2,i] > 0.75
    ]
)
Lar.deleteModelVertices!(model, toDelV)

GL.VIEW([GL.GLGrid(model.G, Lar.cop2lar(model.T[1]))])


#=== NOW IT LOOPS ===#
LarA = Lar.Arrangement

sigma           = spzeros(Int8, 0)
return_edge_map = false
multiproc       = false

model, edge_map = LarA.planar_arrangement_1(model, sigma, true, multiproc)

bicon_comps = LarA.biconnected_components(model.T[1])
V           = convert(Lar.Points, model.G')
copEV       = model.T[1]

# arrangement of isolated components
n = size(bicon_comps, 1)
shells = Array{Lar.Chain, 1}(undef, n)
boundaries = Array{Lar.ChainOp, 1}(undef, n)
EVs = Array{Lar.ChainOp, 1}(undef, n)
# for each component
for p in 2 : n
    ev = copEV[sort(bicon_comps[p]), :]
    # computation of 2-cells
    fe = LarA.minimal_2cycles(V, ev)
    # exterior cycle
    shell_num = LarA.get_external_cycle(V, ev, fe)
    # decompose each fe (co-boundary local to component)
    EVs[p] = ev
    tokeep = setdiff(1:fe.m, shell_num)
    boundaries[p] = fe[tokeep, :]
    shells[p] = fe[shell_num, :]
end

p = 1
ev = copEV[sort(bicon_comps[p]), :]
# computation of 2-cells
fe = LarA.minimal_2cycles(V, ev)                         # <---| Here it loops
# exterior cycle
shell_num = LarA.get_external_cycle(V, ev, fe)
# decompose each fe (co-boundary local to component)
EVs[p] = ev
tokeep = setdiff(1:fe.m, shell_num)
boundaries[p] = fe[tokeep, :]
shells[p] = fe[shell_num, :]
