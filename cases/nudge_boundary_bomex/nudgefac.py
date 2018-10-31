import matplotlib.pyplot as pl
import numpy as np

pl.close('all')

# Numerical grid
xsize = 12800
dx    = 50
x     = np.arange(0.5*dx, xsize, dx)

# Boundary settings
dxb   = 300     # Width of boundary nudging
oxb   = 1500    # Offset from lateral boundary

# Center of nudging area
bc1 = oxb
bc2 = xsize-oxb

f1 = np.exp(-0.5*((x-bc1)/dxb)**2)
f2 = np.exp(-0.5*((x-bc2)/dxb)**2)

pl.figure()
pl.plot(x, f1)
pl.plot(x, f2)
