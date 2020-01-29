# using Lar
using SparseArrays

#=== GEOMETRY ===#
triforce = Lar.Model([
	0.0 4.0 2.0 8.0 6.0 4.0 4.0 3.0 5.0 
	0.0 6.0 3.0 0.0 3.0 0.0 1.0 2.5 2.5 
]);

#=== TOPOLOGY : 1-CELLS ===#
I = Array{Int64}([1, 5, 2, 3, 1, 2, 7, 8, 4, 6, 3, 4, 8, 9, 5, 6, 7, 9, 10, 12, 10, 11, 11, 12, ]);
J = Array{Int64}([1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 9, 9, ]);
K = Array{Int8}([-1, -1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, -1, 1, -1, 1, 1, ]);

Lar.addModelCells!(triforce, 1, SparseArrays.sparse(I, J, K, 12, 9))


#=== TOPOLOGY : 2-CELLS ===#
I = Array{Int64}([2, 3, 3, 4, 2, 4, 1, 2, 1, 3, 1, 4, 1, 5, 1, 5, 1, 5, ]);
J = Array{Int64}([1, 2, 3, 4, 5, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, ]);
K = Array{Int8}([-1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, ]);

Lar.addModelCells!(triforce, 2, SparseArrays.sparse(I, J, K, 5, 12))


