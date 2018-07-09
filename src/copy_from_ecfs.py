import datetime
import os
import subprocess
import tarfile

def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)

def execute(task):
    return_code = subprocess.call(task, shell=True, executable='/bin/bash')
    if return_code != 0:
        raise RuntimeError('subprocess.call failed for {}'.format(task))

if __name__ == '__main__':

    start = datetime.datetime(2016,12,1)
    end   = datetime.datetime(2016,12,10)

    n_cycles = int((end-start).total_seconds() / 3600. / 3.)

    # Path of DDH data:
    path = '/scratch/ms/nl/nkbs/DOWA/LES_forcing'

    # Path in ECFS archive:
    ecfs_path = 'ec:/nkl/harmonie/DOWA/DOWA_40h12tg2_fERA5/ptD_2015-2017/'

    # Copy files from ECFS
    for t in range(n_cycles):
        date = start + t * datetime.timedelta(hours=3)
        print('Processing {}'.format(date))

        # Check if directories exists, if not, create them
        year_path  = '{0:}/{1:04d}'.format(path,       date.year)
        month_path = '{0:}/{1:02d}'.format(year_path,  date.month)
        day_path   = '{0:}/{1:02d}'.format(month_path, date.day)
        cycle_path = '{0:}/{1:02d}'.format(day_path,   date.hour)

        mkdir(year_path)
        mkdir(month_path)
        mkdir(day_path)
        mkdir(cycle_path)

        # Check if DDH tar exists, if not copy from ECFS
        ddh_name = 'DDH_{0:04d}{1:02d}{2:02d}{3:02d}.tar.gz'.format(date.year, date.month, date.day, date.hour)
        f_local  = '{0:}/{1:}'.format(cycle_path, ddh_name)

        ecfs_dir = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}'.format(ecfs_path, date.year, date.month, date.day, date.hour)
        f_remote = '{0:}/{1:}'.format(ecfs_dir, ddh_name)

        if not os.path.exists(f_local):
            print('Copying DDH tar from ECFS')
            execute('ecp {0:} {1:}'.format(f_remote, f_local))

        # Unpack DDH .tar.gz
        if not os.path.exists('{0:}/DHFDLHARM+0010'.format(cycle_path)):
            tar = tarfile.open(f_local)
            tar.extractall(path=cycle_path)
            tar.close()

        # Copy soil temperature/moisture files (daily files)
        if date.hour == 0:
            for patch in ['P01']:
                for level in ['L01', 'L02', 'L03']:
                    f_name = 'wsa_{0:}.{1:}.sfx.NETHERLANDS.DOWA_40h12tg2_fERA5_ptD.{2:04d}{3:02d}{4:02d}.nc'.format(level, patch, date.year, date.month, date.day)
                    f_local = '{0:}/{1:}'.format(cycle_path, f_name)
                    f_remote = '{0:}/{1:}'.format(ecfs_dir, f_name)

                    if not os.path.exists(f_local):
                        print('Copying {}'.format(f_name))
                        execute('ecp {0:} {1:}'.format(f_remote, f_local))

                for level in ['L01', 'L02']:
                    f_name = 'tg_{0:}.{1:}.sfx.NETHERLANDS.DOWA_40h12tg2_fERA5_ptD.{2:04d}{3:02d}{4:02d}.nc'.format(level, patch, date.year, date.month, date.day)
                    f_local = '{0:}/{1:}'.format(cycle_path, f_name)
                    f_remote = '{0:}/{1:}'.format(ecfs_dir, f_name)

                    if not os.path.exists(f_local):
                        print('Copying {}'.format(f_name))
                        execute('ecp {0:} {1:}'.format(f_remote, f_local))
