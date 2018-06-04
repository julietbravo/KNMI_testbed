import numpy as np
import matplotlib.pyplot as pl
from collections import OrderedDict as odict

import datetime
import sys
import os

# Add src directory to Python path
sys.path.append(os.path.abspath('{}/../src/'.format(os.path.dirname(os.path.abspath(__file__)))))

# DALES specific tools to read/write in- and output:
import DALES_tools as dt

pl.close('all')

# Read DALES namelist
nl = dt.Read_namelist('namoptions.001')

# Vertical grid
# ----------------
grid = dt.Linear_stretched_grid(25, nl['domain']['kmax'], 0.02)
#grid = dt.Stretched_grid(nl['domain']['kmax'], 90, 20, dz1=25, dz2=100)
grid.plot()

# Initial profiles
# ----------------
thl = 290   + 6e-3 * grid.z             # Liquid water potential temperature (K)
qt  = 10e-3 - 2e-6 * grid.z             # Total specific humidity (kg m-3)
qt[qt<0] = 0.
u   = np.ones(grid.kmax) * 1             # u-component wind (m s-1)
v   = np.ones(grid.kmax) * 2             # v-component wind (m s-1)
tke = np.ones(grid.kmax) * 0.1           # Subgrid TKE (m2 s-2)

# Write output to prof.inp
# Use ordered dictionary as DALES is input column order sensitive
output = odict({'z (m)':grid.z, 'thl (K)':thl, 'qt (kg kg-1)':qt, 'u (m s-1)':u, 'v (m s-1)':v, 'tke (m2 s-2)':tke})
dt.write_profiles('prof.inp.001', output, 'DOWA Harmonie testbed')

# Initial scalar profiles (for microphysics)
zero = np.zeros(grid.kmax)
output = odict({'z (m)':grid.z, 'qr (kg kg-1)':zero, 'nr (kg kg-1)':zero})
dt.write_profiles('scalar.inp.001', output, 'DOWA Harmonie testbed')

# Time dependent large-scale forcings
# --------------------
t_sfc = np.arange(0,10800.01,300)                   # Time of surface forcing (s)
wthls = np.ones_like(t_sfc)*0.1                     # Surface potential temperature flux (K m s-1)
wqts  = np.ones_like(t_sfc)*0.1e-3                  # Surface specific humidity flux (g kg-1 m s-1)
ps    = np.ones_like(t_sfc)*1e5                     # Surface pressure (pa)
thls  = np.ones_like(t_sfc)*290                     # Surface temperature (K)
qts   = np.ones_like(t_sfc)*10e-3                   # Surface specific humidity (kg kg-1)

t_ls  = np.arange(0,10800.01,3600)                  # Time of atmospheric forcing (s)
ug    = np.ones((t_ls.size, grid.kmax)) * 1           # u-component geostrophic wind (m s-1)
vg    = np.ones((t_ls.size, grid.kmax)) * 2           # v-component geostrophic wind (m s-1)
dtqt  = np.ones((t_ls.size, grid.kmax)) * 1e-3/3600.  # Large-scale specific humidity tendency (kg m-3 s-1)
dtthl = np.ones((t_ls.size, grid.kmax)) * 2/3600.     # Large scale potential temperature tendency (K s-1)
dtu   = np.ones((t_ls.size, grid.kmax)) * 5/3600.     # Large-scale u-component tendency (m s-2)
dtv   = np.ones((t_ls.size, grid.kmax)) * -5/3600.    # Large-scale u-component tendency (m s-2)

# Write to ls_flux.inp
output_sfc = odict({'time':t_sfc, 'wthl_s':wthls, 'wqt_s':wqts, 'p_s':ps, 'thl_s':thls, 'qt_s':qts})
output_ls  = odict({'time':t_ls, 'z':grid.z, 'ug':ug, 'vg':vg, 'dqtdt':dtqt, 'dthldt':dtthl, 'dudt':dtu, 'dvdt':dtv})
dt.write_forcings('ls_flux.inp.001', output_sfc, output_ls, 'DOWA Harmonie testbed')

# Dummy forcings for the microphysics scalars
dt.write_dummy_forcings('ls_fluxsv.inp.001', 2, grid.z)

# Also create non-time dependent input file (lscale.inp)
output_ls  = odict({'height': grid.z, 'ug':ug[0,:], 'vg':vg[0,:], 'wfls':zero, \
                    'dqtdxls':zero, 'dqtdyls':zero, 'dqtdtls': zero, 'dthldt':zero})
dt.write_profiles('lscale.inp.001', output_ls, 'DOWA Harmonie testbed')
