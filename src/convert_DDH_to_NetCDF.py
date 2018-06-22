import datetime
import multiprocessing

from read_Harmonie import *

def convert_cycle(settings):
    # Settings are passed in a dictionary to allow parallel processing
    date = settings['date']
    step = settings['step']
    path = settings['path']

    print('Converting {}'.format(date))

    # Path to the DDH files:
    data_root = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/'.format(path, date.year, date.month, date.day, date.hour)

    # Output NetCDF file:
    out_file = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'.format(path, date.year, date.month, date.day, date.hour)

    # Read DDH files, and convert to NetCDF
    data = Read_DDH_files(data_root, 180, step, quiet=True)
    data.to_netcdf(out_file)


if __name__ == '__main__':
    # Which period to convert?
    start = datetime.datetime(year=2010, month=2, day=28, hour=0)
    end   = datetime.datetime(year=2010, month=2, day=28, hour=21)

    # Path of DDH data. Data structure below is expected to be in format "path/yyyy/mm/dd/hh/"
    path  = '/nobackup/users/stratum/DOWA/LES_forcing'

    # DDH output settings
    step  = 10       # DDH output interval

    # Number of cycles to convert
    n_cycles = int((end-start).total_seconds() / 3600. / 3.) + 1

    # Create list of cycles to convert
    queue = []
    for i in range(n_cycles):
        date = start + i * datetime.timedelta(hours=3)
        queue.append(dict(date=date, step=step, path=path))

    # Timer...
    t0 = datetime.datetime.now()

    # Convert DDH to NetCDF in parallel
    pool = multiprocessing.Pool(processes=4)
    pool.map(convert_cycle, queue)

    dt = datetime.datetime.now() - t0
    print('Conversion of {} files took: {}'.format(len(queue), dt))
