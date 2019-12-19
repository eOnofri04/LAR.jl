using SparseArrays, LightXML
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
