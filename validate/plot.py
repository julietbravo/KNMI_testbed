import matplotlib.pyplot as pl
import datetime
import xarray as xr

pl.close('all')

# Which period to convert?
start = datetime.datetime(year=2016, month=12, day=1, hour=0)
end   = datetime.datetime(year=2016, month=12, day=8, hour=0)

# Path of DDH data. Data structure below is expected to be in format "path/yyyy/mm/dd/hh/"
#path = '/nobackup/users/stratum/DOWA/LES_forcing'
path = '/scratch/ms/nl/nkbs/DOWA/LES_forcing/'

# DDH output settings
step  = 10       # DDH output interval

# Number of cycles to convert
n_cycles = int((end-start).total_seconds() / 3600. / 3.)

# Create list of cycles to convert
files = []
for i in range(n_cycles):
    date = start + i * datetime.timedelta(hours=3)
    files.append('{0}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'\
        .format(path, date.year, date.month, date.day, date.hour))

if 'nc' not in locals():
    nc = xr.open_mfdataset(files)

loc = 2 #7+12
print(nc.name[loc].values[loc])

pl.figure()
ax=pl.subplot(111)
pl.plot(nc.time, nc.H[:,loc], 'k-', label='H')
pl.plot(ax.get_xlim(), [0,0], 'k--')
pl.legend(loc=0)
ax=ax.twinx()
pl.plot(nc.time, nc.T_s[:,loc]-nc['T'][:,loc,0], 'r-', label='Ts-Ta')
pl.plot(ax.get_xlim(), [0,0], 'r--')
pl.legend(loc=2)
