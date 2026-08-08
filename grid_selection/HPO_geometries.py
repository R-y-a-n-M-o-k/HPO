import numpy as np 
N_max = 34 
equilibrium_vals = np.array([1.48332676, 1.45443504, 104.3977314])

#form 1D grids 
r_PO = np.array([0, -0.02, 0.02,
                    -0.03, 0.03, -0.04, 0.04,             
                    -0.05, 0.05, -0.06, 0.07,             
                    -0.07, 0.09, -0.08, 0.11,             
                    -0.10, 0.12, -0.12, 0.15,             
                    -0.14, 0.18, -0.16, 0.21,
                    -0.18, 0.25, -0.20, 0.28,
                    -0.22, 0.33, -0.23, 0.38,
                    -0.24, 0.44, -0.25, 0.50
                ])
r_PO = r_PO + equilibrium_vals[0]
print(f'length of r_PO: {len(r_PO)}')

r_PH = np.array([0, -0.01, 0.01,
                    -0.02, 0.02, -0.03, 0.03,
                    -0.05, 0.05, -0.07, 0.08,
                    -0.09, 0.12, -0.11, 0.16,
                    -0.13, 0.21, -0.16, 0.28,
                    -0.19, 0.37, -0.23, 0.47,
                    -0.27, 0.58, -0.31, 0.72,
                    -0.34, 0.90, -0.37, 1.12,
                    -0.39, 1.36, -0.42, 1.64,
                ])
r_PH = r_PH + equilibrium_vals[1]
print(f'length of r_PH: {len(r_PH)}')

a_HPO = np.array([0, -1,  2,
                     -2,  3,  -4,  5,
                     -6,  8,  -8,  11,
                     -10, 14, -12, 17, 
                     -15, 21, -18, 25,
                     -21, 29, -24, 33,
                     -28, 37, -32, 42,
                     -36, 47, -40, 52, 
                     -44, 57, -47, 62,
        ])
a_HPO = a_HPO + equilibrium_vals[2]
print(f'length of a_HPO: {len(a_HPO)}')

print()
print(r_PO)
print()
print(r_PH)
print()
print(a_HPO)
print()

geometries = []
point = 0

i = 0
while i <= N_max:
    j = 0
    while i + j <= N_max:
        k = 0
        while i + j + k <= N_max:
            geometries.append([point, r_PO[i], r_PH[j], a_HPO[k]])
            point += 1
            k += 1
        j += 1
    i += 1


print('Number of Geometries:')
print(len(geometries))
geometries = np.array(geometries)

np.savetxt('HPO_3D_grid.txt', 
           geometries, 
           fmt=('%d','%.6f','%.6f','%.6f'), 
           delimiter=' ')
