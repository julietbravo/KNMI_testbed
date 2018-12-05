import numpy as np

itot = 2
jtot = 3
ktot = 4

# Create dummy data
a = np.arange(itot*jtot*ktot).reshape((itot,jtot,ktot))

for i in range(itot):
    for j in range(jtot):
        for k in range(ktot):
            print(i,j,k,a[i,j,k])

# Write to binary file
a.T.tofile('bla.bin')
