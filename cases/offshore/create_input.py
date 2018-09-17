import xarray as xr
import numpy as np
import datetime
import sys
import os

from collections import OrderedDict as odict
from scipy.special import erf

# Add src directory to Python path, and import DALES specific tools
src_dir = os.path.abspath('{}/../../src/'.format(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(src_dir)

from DALES_tools import *


# ----- Settings -----
expnr   = 1          # DALES experiment number
iloc    = 0+12       # Location (domain) in NetCDF file
n_accum = 6         # Number of time steps to accumulate in the forcings

# Start and endtime of experiment:
start = datetime.datetime(year=2017, month=6, day=11, hour=0)
end   = datetime.datetime(year=2017, month=6, day=18, hour=0)

path  = '/nobackup/users/stratum/DOWA/LES_forcing'
#path  = '/Users/bart/meteo/data/Harmonie_LES_forcing/'      # Path of DDH data.
#path  = '/scratch/ms/nl/nkbs/DOWA/LES_forcing/'
# ----- End settings -----

# Get list of NetCDF files which need to be processed, and open them with xarray
nc_files = get_file_list(path, start, end)
nc_data  = xr.open_mfdataset(nc_files)

# Get start and end indices in `nc_data`
t0, t1 = get_start_end_indices(start, end, nc_data.time.values)

# Docstring for DALES input files
domain    = str(nc_data.name[0,iloc].values)
lat       = float(nc_data.central_lat[0,iloc].values)
lon       = float(nc_data.central_lon[0,iloc].values)
docstring = '{0} ({1:.2f}N, {2:.2f}E): {3} to {4}'.format(domain, lat, lon, start, end)

print(docstring)

# Create stretched vertical grid for LES
grid = Grid_stretched(kmax=160, dz0=20, nloc1=80, nbuf1=20, dz1=130)
#grid = Grid_stretched(kmax=80, dz0=30, nloc1=40, nbuf1=10, dz1=200)    # debug
#grid = Grid_equidist(kmax=32, dz0=100)  # debug-debug..
#grid.plot()

# Create and write the initial vertical profiles (prof.inp)
create_initial_profiles(nc_data, grid, t0, t1, iloc, docstring, expnr)

# Create and write the surface and atmospheric forcings (ls_flux.inp, ls_fluxsv.inp, lscale.inp)
create_ls_forcings(nc_data, grid, t0, t1, iloc, docstring, n_accum, expnr, harmonie_rad=True)

# Write the nudging profiles (nudge.inp)
z0_nudge = 1500         # Starting height of nudging (m)
dz_nudge = 2000         # Transition thickness
f0_nudge = 0.25         # Nudging near surface
d_nudge  = 1-f0_nudge
nudgefac = (erf((grid.z-z0_nudge)/(0.25*dz_nudge))+1) * d_nudge/2 + f0_nudge  # Nudge factor (0-1)
nudgefac = 1./nudgefac
create_nudging_profiles(nc_data, grid, nudgefac, t0, t1, iloc, docstring, n_accum, expnr)

# Create NetCDF file with reference profiles for RRTMG
create_backrad(nc_data, t0, iloc, expnr)

# Update namelist
namelist = 'namoptions.{0:03d}'.format(expnr)
replace_namelist_value(namelist, 'xlat',  lat)
replace_namelist_value(namelist, 'xlon',  lon)
replace_namelist_value(namelist, 'xday',  start.timetuple().tm_yday)
replace_namelist_value(namelist, 'xtime', start.hour)
replace_namelist_value(namelist, 'kmax',  grid.kmax)
