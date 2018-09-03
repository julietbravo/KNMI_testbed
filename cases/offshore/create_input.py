import numpy as np
import xarray as xr
from collections import OrderedDict as odict
import datetime

# tools.py contains most routines which create/write
# the LES initial conditions, large-scale forcings, ...
from tools import *

# DALES experiment number
expnr = 1

# Location (domain) in NetCDF file
# 0=FINO1, 1=Goeree, 2=EPL, 3=K13, 4=HKZ, 5=P11B, 6=F3FB1
# 7=Cabauw, 8=Loobos, 9=Lutjewad, 10=Schiphol, 11=Rotterdam
# +12 = 10x10km, +24 = 30x30 km
iloc = 0+12

# Start and endtime of experiment:
start = datetime.datetime(year=2017, month=3, day=27, hour=3)
end   = datetime.datetime(year=2017, month=3, day=28, hour=21)

# Path of DDH data. Data structure below is expected to be in format "path/yyyy/mm/dd/hh/"
#path  = '/nobackup/users/stratum/DOWA/LES_forcing'
path  = '/Users/bart/meteo/data/Harmonie_DDH/'

# Get list of NetCDF files which need to be processed, and open them with xarray
nc_files = get_file_list(path, start, end)
nc_data  = xr.open_mfdataset(nc_files)

# Get start and end indices in `nc_data`
t0, t1 = get_start_end_indices(start, end, nc_data.time.values)

# Docstring for DALES input files
domain    = nc_data.name[0,iloc].values
lat       = nc_data.central_lat[0,iloc].values
lon       = nc_data.central_lon[0,iloc].values
docstring = '{0} ({1:.2f}N, {2:.2f}E): {3} to {4}'.format(domain, lat, lon, start, end)

print(docstring)

# Create stretched vertical grid for LES
#grid = Grid_stretched(kmax=160, dz0=20, nloc1=80, nbuf1=20, dz1=130)
grid = Grid_stretched(kmax=80, dz0=30, nloc1=40, nbuf1=10, dz1=200)    # debug
#grid.plot()

# Create and write the initial vertical profiles (prof.inp)
create_initial_profiles(nc_data, grid, t0, t1, iloc, docstring, expnr)

# Create and write the surface and atmospheric forcings (ls_flux.inp, ls_fluxsv.inp, lscale.inp)
create_ls_forcings(nc_data, grid, t0, t1, iloc, docstring, expnr, harmonie_rad=False)

# Write the nudging profiles (nudge.inp)
create_nudging_profiles(nc_data, grid, t0, t1, iloc, docstring, expnr)

# Create NetCDF file with reference profiles for RRTMG
create_backrad(nc_data, t0, iloc, expnr)

# Update namelist
namelist = 'namoptions.{0:03d}'.format(expnr)
replace_namelist_value(namelist, 'xlat',  lat)
replace_namelist_value(namelist, 'xlon',  lon)
replace_namelist_value(namelist, 'xday',  start.timetuple().tm_yday)
replace_namelist_value(namelist, 'xtime', start.hour)
replace_namelist_value(namelist, 'kmax',  grid.kmax)
