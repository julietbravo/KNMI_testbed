import xarray as xr
import datetime
import os

def read_DDH_netcdf(start, end, path):
    """
    Read all converted DDH files in NetCDF format,
    from `start` to `end` date
    """

    print('Reading DDH-NetCDF from {} to {}'.format(start, end))

    nt = int((end-start).total_seconds() / 3600. / 3.)+1

    # Generate list of NetCDF files to read
    files = []
    for t in range(nt):
        date = start + t*datetime.timedelta(hours=3)

        f = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'.\
            format(path, date.year, date.month, date.day, date.hour)

        if os.path.exists(f):
            files.append(f)
        else:
            print('Can not find {}!! Skipping..'.format(f))

    # Read data with xarray
    return xr.open_mfdataset(files)

if __name__ == '__main__':

    path = '/Users/bart/meteo/data/Harmonie_DDH/'

    start = datetime.datetime(year=2017, month=1, day=1,  hour=0)
    end   = datetime.datetime(year=2017, month=1, day=14, hour=21)

    nc = read_DDH_netcdf(start, end, path)
