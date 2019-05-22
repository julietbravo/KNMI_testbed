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

    import argparse

    # Parse input arguments
    p = argparse.ArgumentParser()
    p.add_argument('crosstype', type=str, help='Cross-section type (\"crossxy\", \"crossxz\")')
    p.add_argument('crossname', type=str, help='Variable name (e.g. \"lwp\", \"thl\"), witout \"xy\", \"xz\", ..')
    p.add_argument('expnr',     type=int, help='DALES experiment number')
    p.add_argument('nprocx',    type=int, help='Number of MPI tasks in x-direction')
    p.add_argument('nprocy',    type=int, help='Number of MPI tasks in y-direction')
    p.add_argument('itot',      type=int, help='Number of grid points in x-direction')
    p.add_argument('jtot',      type=int, help='Number of grid points in y-direction')
    p.add_argument('ktot',      type=int, help='Number of grid points in z-direction')
    p.add_argument('--skip',    type=int, help='Only process every \"skip\"-th cross-section in time')
    args = p.parse_args()

    # Some default values (if missing)
    skip = 1 if args.skip is None else args.skip

    # Merge cross-section 
    merge_cross(args.crosstype, args.crossname, args.expnr, args.nprocx, args.nprocy,
                args.itot, args.jtot, args.ktot, skip) 
