
def create_runscript(expname, ntasks, workdir, expnr):
    f = open('run.PBS', 'w')
    f.write('#!/bin/ksh\n')
    f.write('#PBS -S /usr/bin/ksh\n')
    f.write('#PBS -q np\n')
    f.write('#PBS -N {}\n'.format(expname))
    f.write('#PBS -m a\n')
    f.write('#PBS -l walltime=24:10:00\n\n')

    f.write('#PBS -l EC_total_tasks={}\n'.format(ntasks))
    f.write('#PBS -l EC_threads_per_task=1\n')
    f.write('#PBS -l EC_memory_per_task=1GB\n')
    f.write('#PBS -l EC_hyperthreads=1\n\n')

    f.write('prgenvswitchto intel\n')
    f.write('module load netcdf4\n\n')

    # Switch to working directory
    f.write('cd {}\n\n'.format(workdir))

    f.write('aprun -n {0} ./dales4 namoptions.{1:03d} > dales4.out'.format(ntasks, expnr))

    f.close()


def create_postscript(workdir, expnr, itot, jtot, ktot, nprocx, nprocy):
    f = open('post.PBS', 'w')
    f.write('#!/bin/ksh\n')
    f.write('#PBS -S /usr/bin/ksh\n')
    f.write('#PBS -q ns\n')
    f.write('#PBS -N post\n')
    f.write('#PBS -m a\n')
    f.write('#PBS -l walltime=06:00:00\n\n')

    # Switch to working directory
    f.write('cd {}\n\n'.format(workdir))

    settings = '{} {} {} {} {} {}'.format(expnr, nprocx, nprocy, itot, jtot, ktot)

    f.write('python mergecross.py crossxy lwp {}'.format(settings))
    f.write('python mergecross.py crossxy rwp {}'.format(settings))

    f.close()


if __name__ == '__main__':
    create_runscript('asdf', 96, '/home/asdf/asdf', 1)
    create_postscript('/', 1, 192, 192, 160, 12, 8)
