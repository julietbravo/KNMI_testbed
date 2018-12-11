import netCDF4 as nc4
import numpy as np
import os

expnr = 1
base  = 'crossxy'
#base  = 'crossxzspan'
npx   = 2
npy   = 1
itot  = 128
jtot  = 128
ktot  = 32

if 'xz' in base:
    mode = 'xz'
elif 'xy' in base:
    mode = 'xy'

if mode == 'xy':
    chx = int(itot / npx)
    chy = int(jtot / npy)
elif mode == 'xz':
    chx = int(itot / npx)
    chy = 1

for i in range(npx):
    n = 1 if 'span' in base else npy
    for j in range(n):
        print('Processing i={}/{}, j={}/{}'.format(i+1,npx,j+1,npy))

        # Input file name
        fname = '{0:}.x{1:03d}y{2:03d}.{3:03d}.nc'.format(base, i, j, expnr)

        # Slices in x and y dimensions
        sx = np.s_[i*chx:(i+1)*chx]
        sy = np.s_[j*chy:(j+1)*chy]

        if i==0 and j==0:
            # Open first file to copy dimensions et al.
            src = nc4.Dataset(fname)
            dst = nc4.Dataset('{0:}.{1:03d}.nc'.format(base, expnr), 'w')

            # Copy NetCDF attributes and dimensions
            for name in src.ncattrs():
                dst.setncattr(name, src.getncattr(name))
            for name, dim in src.dimensions.items():
                if dim.isunlimited():
                    dst.createDimension(name, None)
                elif name[0] == 'x':
                    dst.createDimension(name, itot)
                elif name[0] == 'y':
                    dst.createDimension(name, jtot)
                elif name[0] == 'z':
                    dst.createDimension(name, ktot)

            # Create variables
            for name, var in src.variables.items():
                dst.createVariable(name, var.datatype, var.dimensions)

            # Copy time
            dst.variables['time'][:]= src.variables['time'][:]

            src.close()

        # Read NetCDF file and write to merged NetCDF
        src = nc4.Dataset(fname)

        # Loop over, and copy data
        for name, var in src.variables.items():
            if name != 'time':
                if name[0] == 'x':
                    dst.variables[name][sx] = src.variables[name][:]
                elif name[0] == 'y':
                    dst.variables[name][sy] = src.variables[name][:]
                elif name[0] == 'z':
                    dst.variables[name][:] = src.variables[name][:]
                else:
                    if mode == 'xy':
                        dst.variables[name][:,sy,sx] = src.variables[name][:]
                    elif mode == 'xz':
                        dst.variables[name][:,:,sx] = src.variables[name][:]

dst.close()
