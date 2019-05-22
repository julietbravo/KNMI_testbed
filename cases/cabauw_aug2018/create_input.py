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
from create_runscript import create_runscript, create_postscript

def execute(task):
    # Execute `task` and return return code
    return subprocess.call(task, shell=True, executable='/bin/bash')

def execute_ret(call):
    # Execute task and return stdout of process
    sp = subprocess.Popen(call, shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
    return sp.stdout.read().decode("utf-8").rstrip('\n')  # jikes!

def fbool(flg):
    # Convert Python bool to Fortran bool
    return '.true.' if flg else '.false.'


if __name__ == '__main__':

    # --------------------
    # Settings
    # --------------------

    expname   = 'cabauw_20160804_20160818_restart'
    expnr     = 1       # DALES experiment number
    iloc      = 7+12    # Location in DDH/NetCDF files (7+12 = 10x10km average Cabauw)
    n_accum   = 1       # Number of time steps to accumulate in the forcings
    warmstart = True    # Run each day/run as a warm start from previous exp

    if expnr == 1:
        # 24 hour runs (cold or warm starts), starting at 00 UTC.
        start  = datetime.datetime(year=2016, month=8, day=4)
        end    = datetime.datetime(year=2016, month=8, day=6)
        #end    = datetime.datetime(year=2016, month=8, day=19)
        dt_exp = datetime.timedelta(hours=24)   # Time interval between experiments
        t_exp  = datetime.timedelta(hours=24)   # Length of experiment
        eps    = datetime.timedelta(hours=1)

    elif expnr == 2:
        # 8 hour runs (cold start), Small domain, high resolution, starting at 17 UTC
        start  = datetime.datetime(year=2016, month=8, day=4, hour=17)
        end    = datetime.datetime(year=2016, month=8, day=19, hour=17)
        dt_exp = datetime.timedelta(hours=24)   # Time interval between experiments
        t_exp  = datetime.timedelta(hours=12)   # Length of experiment
        eps    = datetime.timedelta(hours=1)
    else:
        raise Exception('Undefined experiment number')


    # Paths to the LES forcings, and ERA5/Cabauw for soil initialisation
    host = socket.gethostname()
    if 'cca' in host or 'ccb' in host or 'ecgb' in host:
        # ECMWF CCA/CCB/ECGATE
        path     = '/scratch/ms/nl/nkbs/LES_forcing'
        path_e5  = '/scratch/ms/nl/nkbs/ERA_soil'
        path_out = '/scratch/ms/nl/nkbs/DALES/KNMI_testbed/{}'.format(expname)

    elif 'barts-mbp' in host or 'Barts-MacBook-Pro.local' in host:
        # Macbook
        path     = '/Users/bart/meteo/data/HARMONIE_LES_forcing/'
        path_e5  = '/Users/bart/meteo/data/ERA5/soil/'
        path_out = '/Users/bart/meteo/models/KNMI_testbed/cases/cabauw_aug2018/{}'.format(expname)

    elif 'knmi' in host:
        # KNMI desktop
        path     = '/nobackup/users/stratum/DOWA/LES_forcing/'
        path_e5  = '/nobackup/users/stratum/ERA5/soil/'
        path_out = '/nobackup/users/stratum/KNMI_testbed/cases/cabauw_aug2018/{}'.format(expname)

    else:
        raise Exception('Unknown compute system!')

    # ------------------------
    # End settings
    # ------------------------

    # Create stretched vertical grid for LES
    if expnr == 1:
        # Grid for full diurnal cycle
        #grid = Grid_stretched(kmax=160, dz0=20, nloc1=80, nbuf1=20, dz1=150)
        grid = Grid_stretched(kmax=64, dz0=30, nloc1=40, nbuf1=20, dz1=150)
    elif expnr == 2:
        # High resolution grid for nocturnal runs
        grid = Grid_stretched(kmax=160, dz0=2, nloc1=120, nbuf1=30, dz1=10)
    #grid.plot()

    date = start
    n = 1
    while date < end:
        print('-----------------------')
        print('Starting new experiment')
        print('-----------------------')

        # In case of warm starts, first one is still a cold one..
        start_is_warm = warmstart and n>1
        start_is_cold = not start_is_warm

        # Round start date (first NetCDF file to read) to the 3-hourly HARMONIE cycles
        offset = datetime.timedelta(hours=0) if date.hour%3 == 0 else datetime.timedelta(hours=-date.hour%3)

        # Get list of NetCDF files which need to be processed, and open them with xarray
        nc_files = get_file_list(path, date+offset, date+t_exp+eps)
        nc_data  = xr.open_mfdataset(nc_files)

        # Get start and end indices in `nc_data`
        t0, t1 = get_start_end_indices(date, date + t_exp + eps, nc_data.time.values)

        # Docstring for DALES input files
        domain    = nc_data.name[0,iloc].values
        lat       = float(nc_data.central_lat[0,iloc].values)
        lon       = float(nc_data.central_lon[0,iloc].values)
        docstring = '{0} ({1:.2f}N, {2:.2f}E): {3} to {4}'.format(domain, lat, lon, date, date + t_exp)

        if start_is_cold:
            # Create and write the initial vertical profiles (prof.inp)
            create_initial_profiles(nc_data, grid, t0, t1, iloc, docstring, expnr)

        # Create and write the surface and atmospheric forcings (ls_flux.inp, ls_fluxsv.inp, lscale.inp)
        create_ls_forcings(nc_data, grid, t0, t1, iloc, docstring, n_accum, expnr, harmonie_rad=False)

        # Write the nudging profiles (nudge.inp)
        nudgefac = np.ones_like(grid.z)     # ??
        create_nudging_profiles(nc_data, grid, nudgefac, t0, t1, iloc, docstring, 1, expnr)

        # Create NetCDF file with reference/background profiles for RRTMG
        create_backrad(nc_data, t0, iloc, expnr)

        if start_is_cold:
            # Get the soil temperature and moisture from ERA5
            tsoil   = get_Tsoil_ERA5  (date, 4.9, 51.97, path_e5)
            phisoil = get_phisoil_ERA5(date, 4.9, 51.97, path_e5)

            # Option to re-scale soil moisture content
            soil_in     = soil_med_fine      # ERA5 grid point soil type
            soil_out    = soil_fine          # ~Cabauw soil type
            old_phisoil = phisoil.copy()
            phisoil     = soil_in.rescale(old_phisoil, soil_out)

        # Update namelist
        namelist = 'namoptions.{0:03d}'.format(expnr)
        replace_namelist_value(namelist, 'lwarmstart', fbool(start_is_warm))
        replace_namelist_value(namelist, 'iexpnr',   '{0:03d}'.format(expnr))
        replace_namelist_value(namelist, 'runtime',  t_exp.total_seconds())
        replace_namelist_value(namelist, 'trestart', t_exp.total_seconds())
        replace_namelist_value(namelist, 'xlat',     lat)
        replace_namelist_value(namelist, 'xlat',     lat)
        replace_namelist_value(namelist, 'xlon',     lon)
        replace_namelist_value(namelist, 'xday',     date.timetuple().tm_yday)
        replace_namelist_value(namelist, 'xtime',    date.hour)
        replace_namelist_value(namelist, 'kmax',     grid.kmax)

        if start_is_cold:
            replace_namelist_value(namelist, 'tsoilav',  array_to_string(tsoil))
            replace_namelist_value(namelist, 'phiwav',   array_to_string(phisoil))
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

        # Read back namelist
        nl = Read_namelist('namoptions.{0:03d}'.format(expnr))

        # Copy/move files to work directory
        wdir = '{0}/{1:04d}{2:02d}{3:02d}'.format(path_out, date.year, date.month, date.day)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        # Create SLURM runscript
        print('Creating runscript')
        create_runscript('LES_{}'.format(n), 96, wdir, expnr)

        # Copy/move files to work directory
        exp_str = '{0:03d}'.format(expnr)
        to_copy = ['namoptions.{}'.format(exp_str), 'rrtmg_lw.nc', 'rrtmg_sw.nc', 'dales4',
                   'prof.inp.{}'.format(exp_str), 'scalar.inp.{}'.format(exp_str), 'mergecross.py']
        to_move = ['backrad.inp.{}.nc'.format(exp_str), 'lscale.inp.{}'.format(exp_str),\
                   'ls_flux.inp.{}'.format(exp_str), 'ls_fluxsv.inp.{}'.format(exp_str),\
                   'nudge.inp.{}'.format(exp_str), 'run.PBS', 'post.PBS']

        print('Copying/moving input files')
        for f in to_move:
            shutil.move(f, '{}/{}'.format(wdir, f))
        for f in to_copy:
            shutil.copy(f, '{}/{}'.format(wdir, f))

        if start_is_warm:
            # Link restart files from `prev_wdir` to the current working directory)
            print('Creating symlinks to restart files')

            hh = int(t_exp.total_seconds()/3600)
            mm = int(t_exp.total_seconds()-(hh*3600))

            for i in range(nl['run']['nprocx']):
                for j in range(nl['run']['nprocy']):
                    for ftype in ['d','s','l']:

                        f_in  = '{0}/init{1}{2:03d}h{3:02d}mx{4:03d}y{5:03d}.{6:03d}'\
                                    .format(prev_wdir, ftype, hh, mm, i, j, expnr)
                        f_out = '{0}/init{1}000h00mx{2:03d}y{3:03d}.{4:03d}'\
                                    .format(wdir, ftype, i, j, expnr)

                        if not os.path.islink(f_out):
                            os.symlink(f_in, f_out)

        # Submit task, accounting for job dependencies
        if start_is_warm:
            tid = execute_ret('qsub -W depend=afterok:{} {}/run.PBS'.format(prev_tid, wdir))
            print('Submitted run: {} (depends on: {})'.format(tid, prev_tid))
        else:
            tid = execute_ret('qsub {}/run.PBS'.format(wdir))
            print('Submitted run: {}'.format(tid))

        # Submit post-processing task
        create_postscript(wdir, expnr, nl['domain']['itot'], nl['domain']['jtot'],
                          nl['domain']['kmax'], nl['run']['nprocx'], nl['run']['nprocy'])
        pid = execute_ret('qsub -W depend=afterok:{} {}/post.PBS'.format(tid, wdir))
        print('Submitted postprocessing: {} (depends on: {})'.format(pid, tid))

        # Advance time...
        date += dt_exp
        n += 1
        prev_wdir = wdir
        prev_tid = tid
