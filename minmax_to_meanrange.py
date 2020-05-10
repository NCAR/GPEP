import numpy as np
import os
# additional function: tmin/tmax to tmean/trange
inpath = '/home/gut428/JRA55_day_ds'
outpath = '/home/gut428/JRA55_day_ds'
for y in range(1979,2019):
    print(y)
    infile1 = inpath + '/JRA55_tmin_' + str(y) + '.npz'
    infile2 = inpath + '/JRA55_tmax_' + str(y) + '.npz'
    outfile1 = outpath + '/JRA55_tmean_' + str(y) + '.npz'
    outfile2 = outpath + '/JRA55_trange_' + str(y) + '.npz'
    if os.path.isfile(outfile1) and os.path.isfile(outfile2):
        continue

    tmin = np.load(infile1)
    lattar = tmin['latitude']
    lontar = tmin['longitude']
    tmin = tmin['data']
    tmax = np.load(infile2)
    tmax=tmax['data']
    tmean = (tmin+tmax)/2
    trange = np.abs(tmax-tmin)
    np.savez_compressed(outfile1, data=np.float32(tmean), latitude=lattar, longitude=lontar)
    np.savez_compressed(outfile2, data=np.float32(trange), latitude=lattar, longitude=lontar)