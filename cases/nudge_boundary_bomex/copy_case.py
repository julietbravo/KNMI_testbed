from shutil import copyfile
import sys
import os
import re

def replace_namelist_value(namelist, variable, new_value):
    with open(namelist, 'r') as source:
        lines = source.readlines()
    with open(namelist, 'w') as source:
        for line in lines:
            source.write(re.sub(r'({}).*'.format(variable), r'\1 = {}'.format(new_value), line))


def replace_value(file, old, new):
    with open(file, 'r') as source:
        lines = source.readlines()
    with open(file, 'w') as source:
        for line in lines:
            source.write(re.sub(old, new, line))


def copy(src, dst):
    # Prevent overwriting existing files
    if not os.path.exists(dst):
        copyfile(src, dst)
    else:
        sys.exit('Cannot copy {} to {}: file exists'.format(src, dst))


# Process command line input (old and new expnr)
old = int(sys.argv[1])
new = int(sys.argv[2])

# Copy files
files = ['namoptions', 'prof.inp', 'lscale.inp', 'nudge.inp', 'run.PBS']
for f in files:
    copy('{0:}.{1:03d}'.format(f,old), '{0:}.{1:03d}'.format(f,new))

# Update expnr in namelist and runscript
replace_namelist_value('namoptions.{0:03d}'.format(new), 'iexpnr', new)
replace_value('run.PBS.{0:03d}'.format(new), '{0:03d}'.format(old), '{0:03d}'.format(new))
