import numpy as np
import netCDF4 as nc
from scipy import io
import os, time, sys
from calendar import isleap

# check if netcdf files are corrupted
# path = '/home/gut428/projects/rpp-kshook/gut428/EMDNA_v1/Estimate'
# for y in range(1979,2019):
#     for e in range(1,100):
#         file = '{}/{}/EMDNA_{}.{:03d}.nc4'.format(path, y, y, e)
#         if not os.path.isfile(file):
#             print('File does not exist', y, e)
#         else:
#             try:
#                 d = nc.Dataset(file)
#                 d2 = d['prcp']
#                 d2 = d['tmean']
#                 d2 = d['trange']
#                 d.close()
#             except:
#                 print('File is corrupted', y, e)

yy = int(sys.argv[1])


def cal_mean(data):
    data[data < -100] = np.nan
    ndays = data.shape[0]
    dout = np.zeros([ndays, 2])
    for i in range(ndays):
        di = data[i]
        dout[i, 0] = np.nanmean(di)
        dout[i, 1] = np.sum(~np.isnan(di))
    return dout


# check mean and valid pixels
path = '/home/gut428/projects/rpp-kshook/gut428/EMDNA_v1/Estimate'
outfile = str(yy) + '_check.mat'

if isleap(yy):
    ndays = 366
else:
    ndays = 365

mm = np.zeros([3, ndays, 100])
num = np.zeros([3, ndays, 100])
flag = 0
for y in range(yy, yy + 1):
    for e in range(1, 101):
        print(y, e)
        file = '{}/{}/EMDNA_{}.{:03d}.nc4'.format(path, y, y, e)
        d = nc.Dataset(file)
        d2 = d['prcp'][:].data
        ndays = d2.shape[0]
        dout = cal_mean(d2)
        mm[0, flag:flag + ndays, e-1] = dout[:, 0]
        num[0, flag:flag + ndays, e-1] = dout[:, 1]
        d2 = d['tmean'][:].data
        dout = cal_mean(d2)
        mm[1, flag:flag + ndays, e-1] = dout[:, 0]
        num[1, flag:flag + ndays, e-1] = dout[:, 1]
        d2 = d['trange'][:].data
        dout = cal_mean(d2)
        mm[2, flag:flag + ndays, e-1] = dout[:, 0]
        num[2, flag:flag + ndays, e-1] = dout[:, 1]
    flag = flag + ndays
io.savemat(outfile,{'mm':mm,'num':num})
