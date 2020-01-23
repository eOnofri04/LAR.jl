using SparseArrays, LightXML

function osm2lar(
        filename::String,
        doCrop::Bool = false,
        doFrame::Bool = false
    )::Lar.Model
    xdoc  = parse_file(filename)
    xroot = root(xdoc)

    cropBox   = get_elements_by_tagname(xroot, "bounds")[1]
    node_list = get_elements_by_tagname(xroot, "node"  )
    way_list  = get_elements_by_tagname(xroot, "way"   )

    #=== EXTRACT VERTEX SET AND VERTEX POSITION MAP ===#

    Vidxs = [parse(Int, attribute(v, "id")) for v in node_list]
    Vmax  = max(Vidxs...)
    Vmap  = SparseArrays.spzeros(Int, Vmax)
    Vnum  = length(Vidxs)
    for i = 1 : Vnum   Vmap[Vidxs[i]] = i   end

    #=== EXTRACT (X,Y) AND RESCALE POINTS IN [0, 1] ===#
    lons  = [parse(Float64, attribute(v, "lon")) for v in node_list]
    lmin  = parse(Float64, attribute(cropBox, "minlon"))
    lmax  = parse(Float64, attribute(cropBox, "maxlon"))
    lons  = (lons .- lmin) ./ (lmax - lmin)

    lats  = [parse(Float64, attribute(v, "lat")) for v in node_list]
    lmin  = parse(Float64, attribute(cropBox, "minlat"))
    lmax  = parse(Float64, attribute(cropBox, "maxlat"))
    lats  = (lats .- lmin) ./ (lmax - lmin)

    #=== BUILD LAR MODEL ===#

    model = Lar.Model([lons'; lats'])

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

    if doCrop
        Lar.mapBoxCrop!(model, [0.0 1.0;0.0 1.0], doFrame)
    end
    return model
end

function mapBoxCrop!(
        model::Lar.Model,
        cropBox::Lar.Points,
        doFrame::Bool
    )::Nothing

    toDelV = convert(Array{Int64,1},
        [i  for i = 1 : size(model, 0, 2)
            if model.G[1, i] <= cropBox[1, 1]
            || model.G[1, i] >= cropBox[1, 2]
            || model.G[2, i] <= cropBox[2, 1]
            || model.G[2, i] >= cropBox[2, 2]
        ]
    )
    Lar.deleteModelVertices!(model, toDelV)

    if doFrame
        idx = size(model, 0, 2) + 1
        Lar.addModelVertices!(model, [
            cropBox[1, 1] cropBox[1, 1] cropBox[1, 2] cropBox[1, 2]
            cropBox[2, 1] cropBox[2, 2] cropBox[2, 2] cropBox[2, 1]
        ])
        Lar.addModelCells!(model, 1, [
            [idx, idx+1], [idx+1, idx+2], [idx+2, idx+3], [idx+3, idx]
        ])
    end
    return
end
