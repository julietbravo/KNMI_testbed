import numpy as np
import xarray as xr
from collections import OrderedDict as odict
import datetime
import shutil
import sys
import os

# Add src directory to Python path, and import DALES specific tools
src_dir = os.path.abspath('{}/../../src/'.format(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(src_dir)

from DALES_tools import *
from read_soil_ERA5 import *
from read_soil_Cabauw import *

# --------------------
# Settings
# --------------------
expnr   = 1       # DALES experiment number
iloc    = 7+12    # Location in DDH files (7=Cabauw, 7+12 = 10x10km average Cabauw)
n_accum = 1       # Number of time steps to accumulate in the forcings

start = datetime.datetime(year=2016, month=8, day=4)
end   = datetime.datetime(year=2016, month=8, day=19)
dt    = datetime.timedelta(hours=24)
eps   = datetime.timedelta(hours=1)

# Start and endtime of experiment:
#start = datetime.datetime(year=2016, month=8, day=9, hour=0)
#end   = datetime.datetime(year=2016, month=8, day=10, hour=0)

# Paths to the LES forcings, and ERA5/Cabauw for soil initialisation
#path    = '/scratch/ms/nl/nkbs/DOWA/LES_forcing'
#path_cb = '/scratch/ms/nl/nkbs/DOWA/Cabauw'

path    = '/nobackup/users/stratum/DOWA/LES_forcing'
path_cb = '/nobackup/users/stratum/Cabauw'
path_e5 = '/nobackup/users/stratum/ERA5/soil'

#path     = '/Users/bart/meteo/data/Harmonie_LES_forcing/'
#path_cb  = '/Users/bart/meteo/observations/Cabauw/'
#path_e5  = '/Users/bart/meteo/data/ERA5/soil/'

# ------------------------
# End settings
# ------------------------

date = start
while date < end:

    # Get list of NetCDF files which need to be processed, and open them with xarray
    nc_files = get_file_list(path, date, date+dt+eps)
    nc_data  = xr.open_mfdataset(nc_files)
    
    # Get start and end indices in `nc_data`
    t0, t1 = get_start_end_indices(date, date+dt+eps, nc_data.time.values)
    
    # Docstring for DALES input files
    domain    = nc_data.name[0,iloc].values
    lat       = nc_data.central_lat[0,iloc].values
    lon       = nc_data.central_lon[0,iloc].values
    docstring = '{0} ({1:.2f}N, {2:.2f}E): {3} to {4}'.format(domain, lat, lon, date, date+dt+eps)
    print(docstring)
    
    # Create stretched vertical grid for LES
    #grid = Grid_stretched(kmax=160, dz0=20, nloc1=80, nbuf1=20, dz1=150)
    grid = Grid_stretched(kmax=100, dz0=20, nloc1=40, nbuf1=10, dz1=200)    # debug
    #grid = Grid_stretched(kmax=48,  dz0=20, nloc1=40, nbuf1=10, dz1=200)    # real debug
    #grid.plot()
    
    # Create and write the initial vertical profiles (prof.inp)
    create_initial_profiles(nc_data, grid, t0, t1, iloc, docstring, expnr)
    
    # Create and write the surface and atmospheric forcings (ls_flux.inp, ls_fluxsv.inp, lscale.inp)
    create_ls_forcings(nc_data, grid, t0, t1, iloc, docstring, n_accum, expnr, harmonie_rad=False)
    
    # Write the nudging profiles (nudge.inp)
    nudgefac = np.ones_like(grid.z)
    create_nudging_profiles(nc_data, grid, nudgefac, t0, t1, iloc, docstring, expnr)
    
    # Create NetCDF file with reference profiles for RRTMG
    create_backrad(nc_data, t0, iloc, expnr)
    
    # Get the soil temperature and moisture from ERA5
    tsoil   = get_Tsoil_ERA5  (date, 4.9, 51.97, path_e5)
    phisoil = get_phisoil_ERA5(date, 4.9, 51.97, path_e5)
    
    # Update namelist
    namelist = 'namoptions.{0:03d}'.format(expnr)
    replace_namelist_value(namelist, 'xlat',    lat)
    replace_namelist_value(namelist, 'xlon',    lon)
    replace_namelist_value(namelist, 'xday',    date.timetuple().tm_yday)
    replace_namelist_value(namelist, 'xtime',   date.hour)
    replace_namelist_value(namelist, 'kmax',    grid.kmax)
    replace_namelist_value(namelist, 'tsoilav', array_to_string(tsoil))
    replace_namelist_value(namelist, 'phiwav',  array_to_string(phisoil))
    replace_namelist_value(namelist, 'tsoildeepav', tsoil[-1])  #????

    # Copy/move files to work directory
    wdir = '{0:04d}{1:02d}{2:02d}'.format(date.year, date.month, date.day)
    if not os.path.exists(wdir):
        os.mkdir(wdir)

    to_copy = ['namoptions.001','rrtmg_lw.nc','rrtmg_sw.nc']
    to_move = ['backrad.inp.001.nc','lscale.inp.001','ls_flux.inp.001',\
               'ls_fluxsv.inp.001','nudge.inp.001','prof.inp.001','scalar.inp.001']

    for f in to_move:
        shutil.move(f, '{}/{}'.format(wdir, f))
    for f in to_copy:
        shutil.copy(f, '{}/{}'.format(wdir, f))

    # Advance time...
    date += dt
