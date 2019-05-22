import netCDF4 as nc4
import numpy as np


def merge_cross(base, variable, expnr, npx, npy, itot, jtot, ktot, ntime):

    # Determine some settings
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
    
    # Process cross-sections of each MPI task
    for i in range(npx):
        n = 1 if 'span' in base else npy
        for j in range(n):
            print('Processing i={}/{}, j={}/{}'.format(i+1,npx,j+1,npy))
    
            # Input file name
            fname = '{0:}.x{1:03d}y{2:03d}.{3:03d}.nc'.format(base, i, j, expnr)
    
            # Slices in x and y dimensions
            sx = np.s_[i*chx:(i+1)*chx]
            sy = np.s_[j*chy:(j+1)*chy]
    
            # Create new NetCDF file if first MPI task
            if i==0 and j==0:
                src = nc4.Dataset(fname)
                dst = nc4.Dataset('{0:}.{1:}.{2:03d}.nc'.format(base, variable, expnr), 'w')
    
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
                    if name == 'time' or name[0] == 'x' or name[0] == 'y' or name[0] == 'z' or name.replace(mode, '') == variable:
                        dst.createVariable(name, var.datatype, var.dimensions)
    
                # Copy time
                dst.variables['time'][:]= src.variables['time'][::ntime]
    
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
                        if name.replace(mode, '') == variable:
                            if mode == 'xy':
                                dst.variables[name][:,sy,sx] = src.variables[name][::ntime]
                            elif mode == 'xz':
                                dst.variables[name][:,:,sx] = src.variables[name][::ntime]
    
    dst.close()


if __name__ == '__main__':
    # Only executed is script is directly called

    import sys

    if len(sys.argv) != 10:
        print('Usage: \"python mergecross.py crossname varname expnr nprocx nprocy itot jtot ktot ntime\"')
        print('E.g.:  \"python mergecross.py crossxy lwp 1 12 8 192 192 160 1\"')
    else:    

        base  = str(sys.argv[1])
        var   = str(sys.argv[2])
        exp   = int(sys.argv[3])
        npx   = int(sys.argv[4])
        npy   = int(sys.argv[5])
        itot  = int(sys.argv[6])
        jtot  = int(sys.argv[7])
        ktot  = int(sys.argv[8])
        ntime = int(sys.argv[9])




        merge_cross(base, var, exp, npx, npy, itot, jtot, ktot, ntime) 
