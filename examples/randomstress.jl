
function randomstress(npts::Int = 100, ndgs::Int = 200)

    V = rand(2, npts)
    model = Lar.Model(rand(2, npts))

    cv = unique([rand(1 : npts, 2) for i = 1 : 10*npts])[1 : ndgs*2]
    cv = [[c[1] c[2]] for c in cv if c[1] != c[2]][1 : ndgs]
    I = []
    J = []
    K = ones(Int8, 2*ndgs)
    for i = 1 : ndgs
      push!(J, cv[i]...)
      push!(I, [i; i]...)
    end
    Lar.addModelCells!(model, 1, SparseArrays.sparse(I, J, K, ndgs, npts))

    GL.VIEW([ GL.GLLines(model.G, LAR.cop2lar(model.T[1])) ])


    model1, edge_map = LarA.planar_arrangement_1(model, spzeros(Int8, 0), true)

    W, EW = LAR.fraglines(1.5,1.5,1.5)((model1.G, LAR.cop2lar(model1.T[1])))
    GL.VIEW([ GL.GLLines(W, EW, GL.COLORS[1]) ])

    bicon_comps = LarA.biconnected_components(model1.T[1])

	#V, EVs = LAR.biconnectedComponent((model1.G, LAR.cop2lar(model1.T[1]))) # 2-connected components
	#VV = [[k] for k = 1 : size(V,2)]
	#meshes = [ GL.numbering(.01)((V, (VV,EVs[k])), GL.COLORS[k],1) for k=1:length(EVs) ]
	#GL.VIEW( cat(meshes) );


    edge_map = LarA.remove_mapping_dangling_edges(edge_map, bicon_comps)
    model2 = LarA.planar_arrangement_2(model1, bicon_comps)

    EVs = Lar.FV2EVs(model2)
	FVs = Lar.triangulate2D(model2)
	GL.VIEW(GL.GLExplode(model2.G, FVs, 1.2, 1.2, 1.2, 99, 1));
	GL.VIEW(GL.GLExplode(model2.G, FVs,   1,   1,   1, 99, 1));
	GL.VIEW(GL.GLExplode(model2.G, EVs,   1,   1,   1, 99, 1));

end
