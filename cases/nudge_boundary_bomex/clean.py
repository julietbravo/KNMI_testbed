import glob
import os

all_files = glob.glob('*')
exclude   = ['namoptions.001','namoptions.002','run.PBS']

for f in all_files:
    if f not in exclude and '.py' not in f:
        try:
            os.remove(f)
        except:
            print('Can not remove {}'.format(f))
