import matplotlib.pyplot as pl
import numpy as np

pl.close('all')
pl.ion()

# Numerical grid
xsize = 60000/6.
dx    = 50
x     = np.arange(0.5*dx, xsize, dx)

# Boundary settings
dxb   = 125     # Width of boundary nudging
oxb   = 500    # Offset from lateral boundary

# Center of nudging area
bc1 = oxb
bc2 = xsize-oxb

bc3 = 2*oxb

f1 = np.exp(-0.5*((x-bc1)/dxb)**2)
f2 = np.exp(-0.5*((x-bc2)/dxb)**2)
f3 = np.exp(-0.5*((x-bc3)/dxb)**2)

pl.figure()
pl.plot(x, f1)
pl.plot(x, f2)
pl.plot(x, f3)
