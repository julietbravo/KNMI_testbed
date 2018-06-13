import numpy as np
import matplotlib.pyplot as pl
import xarray as xr
from collections import OrderedDict as odict

import datetime
import sys
import os
import glob

# Add src directory to Python path
sys.path.append(os.path.abspath('{}/../src/'.format(os.path.dirname(os.path.abspath(__file__)))))

# DALES specific tools to read/write in- and output:
import DALES_tools as dt
from DALES_tools import constants

def interp_z(z_Harmonie, z_LES, variable):
    return np.interp(z_LES, z_Harmonie, variable)

def interp_z_time(z_Harmonie, z_LES, variable):
    data = np.zeros((variable.shape[0], z_LES.size))
    for t in range(variable.shape[0]):
        data[t,:] = interp_z(z_Harmonie[t,:], z_LES, variable[t,:])
    return data


if __name__ == '__main__':
    pl.close('all')

    # Read DALES namelist
    nl = dt.Read_namelist('namoptions.001')

    # Create vertical grid
    # ====================
    grid = dt.Stretched_grid(nl['domain']['kmax'], 80, 20, dz1=20, dz2=130)
    grid.plot()

    # Read Harmonie initial conditions & forcings
    files = glob.glob('/nobackup/users/stratum/DOWA/LES_forcing/LES_forcings_20100228*')
    #files = glob.glob('/nobackup/users/stratum/DOWA/LES_forcing/LES_forcings_2010022815*')
    files.sort()
    input = xr.open_mfdataset(files)

    # Create initial profiles
    # =======================
    # Interpolate from Harmonie to LES grid
    p   = interp_z(input['zg'][0,:], grid.z, input['p' ][0,:])
    T   = interp_z(input['zg'][0,:], grid.z, input['T' ][0,:])
    qt  = interp_z(input['zg'][0,:], grid.z, input['q' ][0,:])
    ql  = interp_z(input['zg'][0,:], grid.z, input['ql'][0,:])
    u   = interp_z(input['zg'][0,:], grid.z, input['u' ][0,:])
    v   = interp_z(input['zg'][0,:], grid.z, input['v' ][0,:])
    tke = np.ones(grid.kmax) * 0.1

    # Conversions from Harmonie -> DALES
    exner  = (p / constants['p0'])**(constants['rd']/constants['cp'])
    theta  = T / exner
    thetal = theta - constants['lv'] / (constants['cp'] * exner) * ql

    # Write to prof.inp.001
    output = odict({'z (m)':grid.z, 'thl (K)':thetal, 'qt (kg kg-1)':qt, \
                    'u (m s-1)':u, 'v (m s-1)':v, 'tke (m2 s-2)':tke})
    dt.write_profiles('prof.inp.001', output, 'DOWA Harmonie testbed')

    # Initial scalar profiles (for microphysics) are zero
    zero = np.zeros(grid.kmax)
    output = odict({'z (m)':grid.z, 'qr (kg kg-1)':zero, 'nr (kg kg-1)':zero})
    dt.write_profiles('scalar.inp.001', output, 'DOWA Harmonie testbed')

    # Surface and atmospheric forcings
    # ================================
    ev = 10      # Read every `ev` time steps (Harmonie currently has minute output)
    time_sec = (input.time[::ev] - input.time[0]).values / 1e9

    # Surface variables
    ps   = input['ps' ][::ev].values
    Ts   = input['Tsk'][::ev].values
    qs   = input['qsk'][::ev].values
    dums = np.zeros_like(Ts)    # dummy field with zeros

    # Conversion surface variables
    exner = (ps / constants['p0'])**(constants['rd']/constants['cp'])
    ths = Ts / exner

    # Atmospheric forcings
    dtT  = interp_z_time(input['zg'][::ev], grid.z, input['dtT_dyn' ][::ev,:])
    dtu  = interp_z_time(input['zg'][::ev], grid.z, input['dtu_dyn' ][::ev,:])
    dtv  = interp_z_time(input['zg'][::ev], grid.z, input['dtv_dyn' ][::ev,:])
    dtqv = interp_z_time(input['zg'][::ev], grid.z, input['dtqv_dyn'][::ev,:])
    dtql = interp_z_time(input['zg'][::ev], grid.z, input['dtql_dyn'][::ev,:])
    duma = np.zeros_like(dtT)   # dummy field with zeros

    # Conversions atmospheric forcings
    p_tmp = interp_z_time(input['zg'][::ev], grid.z, input['p'][::ev,:])
    exner = (p_tmp / constants['p0'])**(constants['rd']/constants['cp'])
    # Conversion from dtT->dtth incomplete (pressure contribution missing). Same holds for
    # conversion dtth->dtthl; differentating the theta_l eq. results in another pressure tendency term
    dtth  = dtT / exner
    dtthl = dtth - constants['lv'] / (constants['cp'] * exner) * dtql

    # Write to ls_flux.inp
    output_sfc = odict({'time':time_sec, 'wthl_s':dums, 'wqt_s':dums, 'p_s':ps, 'thl_s':ths, 'qt_s':qs})
    output_ls  = odict({'time':time_sec, 'z':grid.z, 'ug':duma, 'vg':duma, \
                        'dqtdt':dtqv, 'dthldt':dtthl, 'dudt':dtu, 'dvdt':dtv})
    dt.write_forcings('ls_flux.inp.001', output_sfc, output_ls, 'DOWA Harmonie testbed')

    # Dummy forcings for the microphysics scalars
    dt.write_dummy_forcings('ls_fluxsv.inp.001', 2, grid.z)

    # Also create non-time dependent input file (lscale.inp)
    zero = np.zeros_like(grid.z)
    output_ls  = odict({'height': grid.z, 'ug':zero, 'vg':zero, 'wfls':zero, \
                        'dqtdxls':zero, 'dqtdyls':zero, 'dqtdtls': zero, 'dthldt':zero})
    dt.write_profiles('lscale.inp.001', output_ls, 'DOWA Harmonie testbed')
