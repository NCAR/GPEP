# divide yearly files to monthly files to facilitate processing

import numpy as np
from auxiliary_merge import *
import os

y1 = int(sys.argv[1])
y2 = int(sys.argv[2])
print(y1,y2)


path_readowngrid = ['/datastore/GLOBALWATER/CommonData/EMDNA/ERA5_day_ds',  # downscaled gridded data
                   '/datastore/GLOBALWATER/CommonData/EMDNA/MERRA2_day_ds',
                   '/datastore/GLOBALWATER/CommonData/EMDNA/JRA55_day_ds']
outpath = ['/home/gut428/ERA5_day_ds',
           '/home/gut428/MERRA2_day_ds',
           '/home/gut428/JRA55_day_ds']

prefix=['ERA5','MERRA2','JRA55']
vars=['prcp','tmean','trange']

for r in range(len(prefix)):
    for v in range(len(vars)):
        for y in range(y1,y2+1):
            file_yr = path_readowngrid[r] + '/' + prefix[r] + '_ds_' + vars[v] + '_' + str(y) + '.npz'
            if not os.path.isfile(file_yr):
                print('file does not exist')
                continue
            date_list, date_number = m_DateList(y, y, 'ByYear')
            dy=np.load(file_yr)
            dy=dy['data']

            for m in range(12):
                print(prefix[r],vars[v],y,m)
                outfile =  outpath[r] + '/' + prefix[r] + '_ds_' + vars[v] + '_' + str(y*100+m+1) + '.npz'
                if os.path.isfile(outfile):
                    continue
                np.savez_compressed(outfile,data=dy[:,:,date_number['mm']==m+1])