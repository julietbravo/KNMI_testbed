import numpy as np
import xarray as xr
from collections import OrderedDict as odict
import datetime

# tools.py contains most routines which create/write
# the LES initial conditions, large-scale forcings, ...
from tools import *

# Location (domain) in NetCDF file
iloc = 0+12

# Start and endtime of experiment:
start = datetime.datetime(year=2016, month=12, day=1, hour=6)
end   = datetime.datetime(year=2016, month=12, day=1, hour=18)

# Path of DDH data. Data structure below is expected to be in format "path/yyyy/mm/dd/hh/"
#path  = '/nobackup/users/stratum/DOWA/LES_forcing'
path  = '/Users/bart/meteo/data/Harmonie_DDH/'

# Get list of NetCDF files which need to be processed, and open them with xarray
nc_files = get_file_list(path, start, end)
nc_data  = xr.open_mfdataset(nc_files)

# Get start and end indices in `nc_data`
t0, t1 = get_start_end_indices(start, end, nc_data.time.values)

# Docstring for DALES input files
domain    = nc_data.name[0][iloc].values
docstring = '{}: {} to {}'.format(domain, start, end)

# Create stretched vertical grid for LES
grid = Grid_stretched(kmax=160, dz0=20, nloc1=80, nbuf1=20, dz1=130)
#grid.plot()

# Create and write the initial vertical profiles (prof.inp)
create_initial_profiles(nc_data, grid, t0, t1, iloc, docstring)

# Create and write the surface and atmospheric forcings (ls_flux.inp, ls_fluxsv.inp, lscale.inp)
create_ls_forcings(nc_data, grid, t0, t1, iloc, docstring)

# Write the nudging profiles (nudge.inp)
create_nudging_profiles(nc_data, grid, t0, t1, iloc, docstring)
