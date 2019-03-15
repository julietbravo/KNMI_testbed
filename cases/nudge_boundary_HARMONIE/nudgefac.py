import matplotlib.pyplot as pl
import numpy as np

pl.close('all')

# Numerical grid
xsize = 64*500
ysize = 64*500
dx    = 200
dy    = 200
x     = np.arange(0.5*dx, xsize, dx)
y     = np.arange(0.5*dy, ysize, dy)

# Boundary settings
dxb   = 800     # Width of boundary nudging
oxb   = 3000    # Offset from lateral boundary
rb    = 2300    # Radius or corners domain

# Center of nudging area
xbc1 = oxb
xbc2 = xsize-oxb

ybc1 = oxb
ybc2 = ysize-oxb

dc = rb+oxb

f = np.zeros((y.size, x.size), dtype=np.float)

def corner_factor(x, y, xc, yc, r, dxb):
    D = np.sqrt((x-xc)**2 + (y-yc)**2) -rb
    return np.exp(-0.5*(D/dxb)**2)

def normal_factor(x, y, xc, yc, r, dxb):
    D = np.sqrt((x-xc)**2 + (y-yc)**2) -rb
    return np.exp(-0.5*(D/dxb)**2)

for i in range(x.size):
    for j in range(y.size):

        if y[j] < dc and x[i] < dc:
            # SW-corner
            f[j,i] += corner_factor(x[i], y[j], dc, dc, rb, dxb)

        elif y[j] < dc and x[i] > xsize-dc:
            # SE-corner
            f[j,i] += corner_factor(x[i], y[j], xsize-dc, dc, rb, dxb)

        elif y[j] > ysize-dc and x[i] < dc:
            # NW-corner
            f[j,i] += corner_factor(x[i], y[j], dc, ysize-dc, rb, dxb)

        elif y[j] > ysize-dc and x[i] > xsize-dc:
            # NE-corner
            f[j,i] += corner_factor(x[i], y[j], xsize-dc, ysize-dc, rb, dxb)

        else:
            f[j,i] += np.exp(-0.5*((x[i]-xbc1)/dxb)**2)
            f[j,i] += np.exp(-0.5*((x[i]-xbc2)/dxb)**2)
            f[j,i] += np.exp(-0.5*((y[j]-ybc1)/dxb)**2)
            f[j,i] += np.exp(-0.5*((y[j]-ybc2)/dxb)**2)


pl.figure()
ax=pl.subplot(121, aspect='equal')
pl.pcolormesh(x,y,f)
pl.colorbar()

ax=pl.subplot(122)
pl.plot(x,f[int(y.size/2),:])

#
#bc3 = 2*oxb
#
#f1 = np.exp(-0.5*((x-bc1)/dxb)**2)
#f2 = np.exp(-0.5*((x-bc2)/dxb)**2)
#f3 = np.exp(-0.5*((x-bc3)/dxb)**2)
#
#pl.figure()
#pl.plot(x, f1)
#pl.plot(x, f2)
#pl.plot(x, f3)
