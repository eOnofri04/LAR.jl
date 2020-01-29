## FROM CLOUD POINT

using SparseArrays
@using LinearAlgebraicRepresentation as LAR
@using ViewerGL as GL

include("./../src/Lar.jl")
LarA = Lar.Arrangement

include("./models/casa.jl")


Lar.addModelCells!(casa, 2, LarA.biconnected_components(casa.T[1]))

bbox = LarA.getModelBoundingBoxes(casa, 2)
bboxXL = [box .+ [-0.5 +0.5] for box in bbox]

model = Lar.Model(zeros(3, 0));

for comp = 1 : length(bboxXL)
    pts = Lar.getModelCellVertices(casa, 2, comp);
    # plane : ax + by + cz = d
    d = 1
    a, b, c = [pts[1]'; pts[floor(Int, 1+end/2)]'; pts[end]']\ones(3, 1)*d

    x = bboxXL[comp][1, :]
    y = bboxXL[comp][2, :]
    z = bboxXL[comp][3, :]

    CM = Lar.Model(zeros(3, 0));
    V = [];

    px = zeros(2, 2)
    px[1, 1] = (d - b*y[1] - c*z[1]) / a
    px[1, 2] = (d - b*y[1] - c*z[2]) / a
    px[2, 1] = (d - b*y[2] - c*z[1]) / a
    px[2, 2] = (d - b*y[2] - c*z[2]) / a
    for j = 1 : 2  for i = 1 : 2  if x[1] < px[i, j] < x[2]
            push!(V, [px[i, j]; y[i]; z[j]])
            Lar.addModelVertex!(CM, [px[i, j]; y[i]; z[j]])
    end  end  end

    py = zeros(2, 2)
    py[1, 1] = (d - a*x[1] - c*z[1]) / b
    py[1, 2] = (d - a*x[1] - c*z[2]) / b
    py[2, 1] = (d - a*x[2] - c*z[1]) / b
    py[2, 2] = (d - a*x[2] - c*z[2]) / b
    for j = 1 : 2  for i = 1 : 2  if y[1] < py[i, j] < y[2]
            push!(V, [x[i]; py[i, j]; z[j]])
            Lar.addModelVertex!(CM, [x[i]; py[i, j]; z[j]])
    end  end  end

    pz = zeros(2, 2)
    pz[1, 1] = (d - a*x[1] - b*y[1]) / c
    pz[1, 2] = (d - a*x[1] - b*y[2]) / c
    pz[2, 1] = (d - a*x[2] - b*y[1]) / c
    pz[2, 2] = (d - a*x[2] - b*y[2]) / c
    for j = 1 : 2  for i = 1 : 2  if z[1] < pz[i, j] < z[2]
        push!(V, [x[i]; y[j]; pz[i, j]])
        Lar.addModelVertex!(CM, [x[i]; y[j]; pz[i, j]])
    end  end  end

    for i = 1 : size(CM, 0, 2)   for j = i + 1 : size(CM, 0, 2)
        if sum(CM.G[:, i] .== CM.G[:, j]) >= 1
            Lar.addModelCell!(CM, 1, [i, j])
        end
    end  end
    Lar.addModelCell!(CM, 2, collect(1:size(CM, 0, 2)))

    Lar.uniteModels!(model, CM)
end
