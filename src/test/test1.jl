Lara = Lar.Arrangement;
V = [
    1.0 0.0 0.0 0.5 1.0;
    0.0 1.0 0.5 1.0 1.0
];
model = Lar.Model(V);
EV = [[1, 2], [2, 5], [3, 4], [4, 5]];
model.T[1] = Lar.coboundary_0(EV::Lar.Cells);
bigPI = Lar.spaceindex((V, EV));
Lara.frag_edge(model, 1, bigPI)[1]
