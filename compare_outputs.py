import numpy as np
import netCDF4 as nc
import auxiliary as au
import regression as reg
import datetime as dt
from matplotlib import pyplot as plt
from scipy import io
import os
import sys
from scipy.interpolate import griddata

file = '/Users/localuser/GMET/pyGMET_res2/regress_daily_output.npz'
dd = np.load(file)
pcp1 = dd['pcp']

file = '/Users/localuser/GMET/Example_tgq/outputs_eCAImodify/regress_daily.nc'
dd = nc.Dataset(file)
pcp2 = dd['pcp'][:].data
dd.close()
pcp2 = np.transpose(pcp2, [1, 2, 0])
pcp2[pcp2 == 0] = -3

pcp1m0 = np.zeros(31)
pcp2m0 = np.zeros(31)
for i in range(31):
    pcp1m0[i] = np.mean(pcp1[:, :, i])
    pcp2m0[i] = np.mean(pcp2[:, :, i])

pcp1 = au.retransform(pcp1, 4, 'box-cox')
pcp2 = au.retransform(pcp2, 4, 'box-cox')

pcp1m = np.zeros(31)
pcp2m = np.zeros(31)
for i in range(31):
    pcp1m[i] = np.mean(pcp1[:, :, i])
    pcp2m[i] = np.mean(pcp2[:, :, i])


# read original station data
path='/Users/localuser/GMET/Example_tgq/StnDaily_train'
files = os.listdir(path)
pcpstn = np.nan * np.zeros([len(files),31])
for i in range(len(files)):
    filef = path + '/' + files[i]
    fid = nc.Dataset(filef)
    if 'prcp' in fid.variables.keys():
        pcpstn[i,:] = fid['prcp'][:].data
    fid.close()
pcpstnm = np.nanmean(pcpstn,axis=0)

plt.figure(figsize=(10, 10))

plt.subplot(221)
plt.imshow(pcp1[:, :, 21])
plt.colorbar()
plt.title('PyGMET Day 22')
plt.clim([0, 5])
plt.subplot(222)
plt.imshow(pcp2[:, :, 21])
plt.title('eCAI Day 22')
plt.clim([0, 5])
plt.colorbar()
plt.subplot(223)
plt.plot(pcp1m0)
plt.plot(pcp2m0)
plt.legend(['PyGMET', 'eCAI'])
plt.xlabel('Days')
plt.ylabel('pcp')
plt.subplot(224)
plt.plot(pcp1m)
plt.plot(pcp2m)
plt.plot(pcpstnm)
plt.legend(['PyGMET', 'eCAI','Station'])
plt.xlabel('Days')
plt.ylabel('pcp (mm)')
plt.show()