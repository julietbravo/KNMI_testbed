import glob
import datetime
import os
import shutil
import tarfile

def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)

data_root = '/nobackup/users/stratum/DOWA/LES_forcing/'

files = glob.glob('{}/tars/*.tar.gz'.format(data_root))
files.sort()

for file in files:
    file_name = file.split('/')[-1]                 # File name minus path
    dtg = file_name.split('_')[-1].split('.')[0]    # yyyymmddhh part

    print('Unpacking {}'.format(file_name))

    date = datetime.datetime.strptime(dtg, '%Y%m%d%H')

    # Absolute path to output directories
    year_dir  = '{0:}/{1:04d}'.format(data_root, date.year)
    month_dir = '{0:}/{1:04d}/{2:02d}'.format(data_root, date.year, date.month)
    day_dir   = '{0:}/{1:04d}/{2:02d}/{3:02d}'.format(data_root, date.year, date.month, date.day)
    hour_dir  = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/'.format(data_root, date.year, date.month, date.day, date.hour)

    # Make sure output directory exists
    mkdir(year_dir )
    mkdir(month_dir)
    mkdir(day_dir  )
    mkdir(hour_dir )

    # Move tar.gz file
    shutil.copy2(file, hour_dir)

    # Unpack
    tar = tarfile.open(file)
    tar.extractall(path=hour_dir)
    tar.close()
