import numpy as np
import auxiliary as au
import netCDF4 as nc
import os
from matplotlib import pyplot as plt
from scipy import io
import seaborn as sns

# 1. read and save
# file = '/Users/localuser/GMET/JRA55_prcp_2018.nc4'
# outfile = '/Users/localuser/GMET/JRA55_prcp_2018_sub.npz'
#
# ncfid = nc.Dataset(file)
# pcp = ncfid['data'][:].data
# pcp = np.transpose(pcp,[2,1,0])
#
# latrange = [5,85]
# lonrange = [-180,-50]
#
# lattar = [35, 40]
# lontar = [-110, -105]
# nrows = int((lattar[1]-lattar[0])/0.1)
# ncols = int((lontar[1]-lontar[0])/0.1)
#
# rowsta = int((latrange[1]-lattar[1])/0.1)
# colsta = int((lontar[0]-lonrange[0])/0.1)
#
# pcptar = pcp[rowsta:rowsta+nrows, colsta:colsta+ncols, :]
#
# np.savez_compressed(outfile, pcprea = pcptar)

# 2. compare the CDF of reanalysis and station data
# read reanalysis data
file1 = '/Users/localuser/GMET/ERA5_prcp_2018_sub.npz'
file2 = '/Users/localuser/GMET/MERRA2_prcp_2018_sub.npz'
file3 = '/Users/localuser/GMET/JRA55_prcp_2018_sub.npz'
pcpera = np.load(file1)
pcpera = pcpera['pcprea']
pcpmerra = np.load(file2)
pcpmerra = pcpmerra['pcprea']
pcpjra = np.load(file3)
pcpjra = pcpjra['pcprea']
pcpera = np.flipud(pcpera)
pcpmerra = np.flipud(pcpmerra)
pcpjra = np.flipud(pcpjra)

# read station info and grid info
FileStnInfo = '/Users/localuser/GMET/pyGMET_exp/inputs/stnlist_example.txt'
FileGridInfo = '/Users/localuser/GMET/pyGMET_exp/inputs/gridinfo_example.nc'
stnID, stninfo = au.readstnlist(FileStnInfo)

ncfid = nc.Dataset(FileGridInfo)
gridlat = ncfid.variables['latitude'][:].data
gridlon = ncfid.variables['longitude'][:].data
gridele = ncfid.variables['elev'][:].data
gridgns = ncfid.variables['gradient_n_s'][:].data
gridgwe = ncfid.variables['gradient_w_e'][:].data
mask = ncfid.variables['mask'][:].data  # 1: grids to be considered; the other values: invalid grids
ncfid.close()

# read station data
FileStnData = '/Users/localuser/GMET/pyGMET_exp/station_data.npz'
datatemp = np.load(FileStnData)
prcp_stn_daily = datatemp['prcp_stn_daily']


nrows, ncols = np.shape(gridlat)
gridinfo = np.zeros([nrows, ncols, 6])
gridinfo[:, :, 0] = 1
gridinfo[:, :, 1] = gridlat
gridinfo[:, :, 2] = gridlon
gridinfo[:, :, 3] = gridele
gridinfo[:, :, 4] = gridgns
gridinfo[:, :, 5] = gridgwe
del gridlat, gridlon, gridele, gridgns, gridgwe

# extract reanalysis data corresponding to stations
nstn = len(stnID)
gridlat = gridinfo[:, 1, 1]
gridlon = gridinfo[1, :, 2]

pcpera2 = np.zeros([nstn,365])
pcpmerra2 = np.zeros([nstn,365])
pcpjra2 = np.zeros([nstn,365])
for i in range(nstn):
    stnlat = stninfo[i, 1]
    stnlon = stninfo[i, 2]
    row = np.argmin(np.abs(stnlat - gridlat))
    col = np.argmin(np.abs(stnlon - gridlon))

    pcpera2[i,:] = pcpera[row, col,:]
    pcpmerra2[i,:] = pcpmerra[row, col, :]
    pcpjra2[i,:] = pcpjra[row, col, :]




kgeera = np.zeros([nstn,4])
kgemerra = np.zeros([nstn,4])
kgejra = np.zeros([nstn,4])
for i in range(nstn):
    kgeera[i,:] = au.kge2012(prcp_stn_daily[i,:], pcpera2[i,:])
    kgemerra[i, :] = au.kge2012(prcp_stn_daily[i, :], pcpmerra2[i, :])
    kgejra[i, :] = au.kge2012(prcp_stn_daily[i, :], pcpjra2[i, :])

# simple merge
w1 = kgeera[:,0].copy()
w2 = kgemerra[:,0].copy()
w3 = kgejra[:,0].copy()
w1[w1<0] = 0
w2[w2<0] = 0
w3[w3<0] = 0
w1 = w1 ** 2
w2 = w2 ** 2
w3 = w3 ** 2
wsum = w1 + w2 + w3
w1 = w1 / wsum
w2 = w2 / wsum
w3 = w3 / wsum
w1[wsum == 0] = 1/3
w2[wsum == 0] = 1/3
w3[wsum == 0] = 1/3

pcpmerge = np.zeros([nstn,365])
for i in range(nstn):
    pcpmerge[i,:] = w1[i]*pcpera2[i,:] + w2[i]*pcpmerra2[i,:]+w3[i]*pcpjra2[i,:]

kgemerge = np.zeros([nstn,4])
for i in range(nstn):
    kgemerge[i, :] = au.kge2012(prcp_stn_daily[i, :], pcpmerge[i, :])

print(np.nanmedian(kgeera,axis=0))
print(np.nanmedian(kgemerra,axis=0))
print(np.nanmedian(kgejra,axis=0))
print(np.nanmedian(kgemerge,axis=0))