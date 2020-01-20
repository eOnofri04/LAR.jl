using LinearAlgebraicRepresentation
LAR = LinearAlgebraicRepresentation

function triangulate2D(model::Lar.Model)::Array{Lar.Cells}
    return LAR.triangulate2D(convert(Lar.Points, model.G'), [model.T[1], model.T[2]])
end

function FV2EVs(model::Lar.Model)::Array{Lar.Cells}
    return LAR.FV2EVs(model.T[1], model.T[2])
end
