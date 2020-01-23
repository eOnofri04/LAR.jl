include("./../src/Lar.jl")
using SparseArrays, ViewerGL, LinearAlgebraicRepresentation
LAR = LinearAlgebraicRepresentation
LarA = Lar.Arrangement
GL = ViewerGL

#=== LAURENTINA DATASET ===#
#  pre-crop E, V = 24043×22808
#  cropped  E, V = 21454×20389
#  arranged E, V =
filename = "examples/osm/laurentina.osm"
filename = "examples/osm/laurentina2.osm"

#=== MARCONI DATASET ===#
#  pre-crop V =     ; E =     ; F =   0
#  cropped  V = 5153; E = 5489; F =   0
#  arranged V = 4600; E = 5062; F = 838
#filename = "examples/osm/marconi.osm"

#=== ALBANO LAZIALE DATASET ===#
#  pre-crop E, V =
#  cropped  E, V = 9205×8879
#  arranged E, V = 9305×8892
#filename = "examples/osm/albano.osm"

model = Lar.osm2lar(filename, true, false)

GL.VIEW([GL.GLGrid(model.G, LAR.cop2lar(model.T[1]), GL.COLORS[3])])

model, edge_map = LarA.planar_arrangement(model, spzeros(Int8, 0), true, false)

GL.VIEW([GL.GLGrid(model.G, LAR.cop2lar(model.T[1]), GL.COLORS[3])])

FVs = Lar.triangulate2D(model)
GL.VIEW(GL.GLExplode(model.G,FVs, 1., 1., 1.,99,1));
GL.VIEW(GL.GLExplode(model.G,FVs,1.2,1.2,1.2,99,1));

GL.VIEW( GL.GLExplode(
    model.G,
    [f for f in FVs if length(f) > 20],
    1.5, 1.5, 1.5, 99, 1)
);

#==

#=== NOW IT LOOPS ===#
LarA = Lar.Arrangement

sigma           = SparseArrays.spzeros(Int8, 0)
return_edge_map = false
multiproc       = false

# model, edge_map = LarA.planar_arrangement_1(model, sigma, true, multiproc)

edgenum = size(model, 1, 1)
edge_map = Array{Array{Int, 1}, 1}(undef,edgenum)
rV = Lar.Points(zeros(0, 2))
rEV = SparseArrays.spzeros(Int8, 0, 0)
finalcells_num = 0

# spaceindex computation
bigPI = LarA.spaceIndex(model)
...

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


==#
