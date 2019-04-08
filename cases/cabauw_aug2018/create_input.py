import numpy as np
import xarray as xr
from collections import OrderedDict as odict
import datetime
import shutil
import sys
import os
import subprocess
import socket

# Add src directory to Python path, and import DALES specific tools
src_dir = os.path.abspath('{}/../../src/'.format(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(src_dir)

from DALES_tools import *
from IFS_soil import *

def execute(task):
    subprocess.call(task, shell=True, executable='/bin/bash')

if __name__ == '__main__':

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

    # Paths to the LES forcings, and ERA5/Cabauw for soil initialisation
    host = socket.gethostname()
    if 'cca' in host or 'ccb' in host or 'ecgb' in host:
        # ECMWF CCA/CCB/ECGATE
        path     = '/scratch/ms/nl/nkbs/LES_forcing'	# CCA/CCB
        path_e5  = '/scratch/ms/nl/nkbs/ERA_soil'	        # CCA/CCB
        path_out = '/scratch/ms/nl/nkbs/DALES/KNMI_testbed/cabauw_20160804_20160818_soil_scaled/'

    elif 'barts-mbp' in host:
        # Macbook
        path     = '/Users/bart/meteo/data/Harmonie_LES_forcing/'	# Macbook
        path_e5  = '/Users/bart/meteo/data/ERA5/soil/'		# Macbook
        path_out = '/Users/bart/meteo/models/KNMI_testbed/cases/cabauw_aug2018/'

    # ------------------------
    # End settings
    # ------------------------

    date = start
    n = 1
    while date < end:

        # Get list of NetCDF files which need to be processed, and open them with xarray
        nc_files = get_file_list(path, date, date+dt+eps)
        nc_data  = xr.open_mfdataset(nc_files)

        # Get start and end indices in `nc_data`
        t0, t1 = get_start_end_indices(date, date+dt+eps, nc_data.time.values)

        # Docstring for DALES input files
        domain    = nc_data.name[0,iloc].values
        lat       = float(nc_data.central_lat[0,iloc].values)
        lon       = float(nc_data.central_lon[0,iloc].values)
        docstring = '{0} ({1:.2f}N, {2:.2f}E): {3} to {4}'.format(domain, lat, lon, date, date+dt+eps)
        print(docstring)

        # Create stretched vertical grid for LES
        grid = Grid_stretched(kmax=160, dz0=20, nloc1=80, nbuf1=20, dz1=150)
        #grid = Grid_stretched(kmax=100, dz0=20, nloc1=40, nbuf1=10, dz1=200)    # debug
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

        # Option to re-scale soil moisture content
        soil_in  = soil_med_fine      # ERA5 grid point soil type
        soil_out = soil_fine          # ~Cabauw soil type
        old_phisoil = phisoil.copy()
        phisoil = soil_in.rescale(old_phisoil, soil_out)

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

        print('Setting soil properties for {} (input={})'.format(soil_out.name, soil_in.name))
        replace_namelist_value(namelist, 'gammasat', soil_out.gammasat)
        replace_namelist_value(namelist, 'nvg',      soil_out.nvg)
        replace_namelist_value(namelist, 'Lvg',      soil_out.lvg)
        replace_namelist_value(namelist, 'alphavg',  soil_out.alphavg)
        replace_namelist_value(namelist, 'phir',     soil_out.phir)
        replace_namelist_value(namelist, 'phi',      soil_out.phi_sat)
        replace_namelist_value(namelist, 'phiwp',    soil_out.phi_wp)
        replace_namelist_value(namelist, 'phifc',    soil_out.phi_fc)

        # Copy/move files to work directory
        wdir = '{0}/{1:04d}{2:02d}{3:02d}'.format(path_out, date.year, date.month, date.day)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        # AARGH, change directory in run script (why?)
        with open('run.PBS') as f:
            lines = f.readlines()
        with open('run.PBS', 'w') as f:
            for line in lines:
                if line[:3] == 'cd ':
                    f.write('cd {}\n'.format(wdir))
                elif '-N' in line:
                    f.write('#PBS -N LES_CB_{}\n'.format(n))
                else:
                    f.write(line)

        to_copy = ['namoptions.001','rrtmg_lw.nc','rrtmg_sw.nc','dales4','run.PBS']
        to_move = ['backrad.inp.001.nc','lscale.inp.001','ls_flux.inp.001',\
                   'ls_fluxsv.inp.001','nudge.inp.001','prof.inp.001','scalar.inp.001']

        for f in to_move:
            shutil.move(f, '{}/{}'.format(wdir, f))
        for f in to_copy:
            shutil.copy(f, '{}/{}'.format(wdir, f))

        # Submit task!
        execute('qsub {}/run.PBS'.format(wdir))

        # Advance time...
        date += dt
        n += 1
