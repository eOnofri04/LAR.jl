# using Lar
using SparseArrays

#=== GEOMETRY ===#
cubes = Lar.Model([
	-0.0797872 -0.1833684 0.0219132 -0.0816681 0.4465752 0.3429939 0.5482755 0.4446942 0.7187785 0.5249728 0.5799289 0.3861232 0.8613671 0.6675614 0.7225174 0.5287117 -0.2241867 0.1406597 -0.0449854 0.319861 -0.12314 0.2417064 0.0560613 0.4209077 
	0.2943071 0.3960075 0.8243972 0.9260975 0.2118999 0.3136002 0.7419899 0.8436902 0.5392904 0.4004408 0.7762076 0.6373579 0.5812714 0.4424218 0.8181886 0.6793389 0.9338056 1.1130069 0.7580241 0.9372254 0.5985085 0.7777098 0.422727 0.6019283 
	-0.1823453 0.344017 -0.2647526 0.2616098 -0.0628419 0.4635205 -0.1452491 0.3811132 0.1230172 0.2656057 0.1649982 0.3075868 0.3577031 0.5002916 0.3996841 0.5422726 0.4075841 0.5086308 0.0722869 0.1733336 0.6373704 0.7384171 0.3020733 0.40312 
]);

#=== TOPOLOGY : 1-CELLS ===#
I = Array{Int64}([1, 5, 9, 1, 6, 10, 2, 5, 11, 2, 6, 12, 3, 7, 9, 3, 8, 10, 4, 7, 11, 4, 8, 12, 13, 17, 21, 13, 18, 22, 14, 17, 23, 14, 18, 24, 15, 19, 21, 15, 20, 22, 16, 19, 23, 16, 20, 24, 25, 29, 33, 25, 30, 34, 26, 29, 35, 26, 30, 36, 27, 31, 33, 27, 32, 34, 28, 31, 35, 28, 32, 36, ]);
J = Array{Int64}([1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, ]);
K = Array{Int8}([-1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, ]);

Lar.addModelCells!(cubes, 1, SparseArrays.sparse(I, J, K, 36, 24))


#=== TOPOLOGY : 2-CELLS ===#
I = Array{Int64}([1, 3, 1, 4, 2, 3, 2, 4, 1, 5, 1, 6, 2, 5, 2, 6, 3, 5, 3, 6, 4, 5, 4, 6, 7, 9, 7, 10, 8, 9, 8, 10, 7, 11, 7, 12, 8, 11, 8, 12, 9, 11, 9, 12, 10, 11, 10, 12, 13, 15, 13, 16, 14, 15, 14, 16, 13, 17, 13, 18, 14, 17, 14, 18, 15, 17, 15, 18, 16, 17, 16, 18, ]);
J = Array{Int64}([1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35, 36, 36, ]);
K = Array{Int8}([1, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, ]);

Lar.addModelCells!(cubes, 2, SparseArrays.sparse(I, J, K, 18, 36))


#=== TOPOLOGY : 3-CELLS ===#
I = Array{Int64}([]);
J = Array{Int64}([]);
K = Array{Int8}([]);

Lar.addModelCells!(cubes, 3, SparseArrays.sparse(I, J, K, 0, 18))


