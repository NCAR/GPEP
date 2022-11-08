# Calcaulate the ratio of undercatch correction (actually this is an overall correction that does not distinguish the type of error)
# this script contains a module to correct undercatch bias using climatology from http://www.gloh2o.org/pbcor/

import numpy as np
import netCDF4 as nc
import sys, os


# ### Mac settings
# path_oi = '/Users/localuser/Research/EMDNA/oimerge'
# near_file_GMET = '/Users/localuser/Research/EMDNA/regression/weight_nearstn.npz'
# gmet_stndatafile = '/Users/localuser/Research/EMDNA/stndata_whole.npz'
# FileStnInfo = '/Users/localuser/Research/EMDNA/basicinfo/stnlist_whole.txt'  # station basic information (lists)
# FileGridInfo = '/Users/localuser/Research/EMDNA/basicinfo/gridinfo_whole.nc'  # study area information
# FileClimo= '/Users/localuser/Research/PBCOR_V1.0/WorldClim_V2.nc'
# ### Mac settings
# FileClimo_tar = '/Users/localuser/Research/PBCOR_V1.0/WorldClim_V2_NA.npz'

### Plato settings
path_oi = '/datastore/GLOBALWATER/CommonData/EMDNA_new/OImerge_GWRLSBMA'
FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/gridinfo_whole.nc'  # study area information
FileClimo= '/datastore/GLOBALWATER/CommonData/EMDNA_new/PBCOR_V1.0/WorldClim_V2.nc'
### Plato settings

FileClimo_tar = '/datastore/GLOBALWATER/CommonData/EMDNA_new/PBCOR_V1.0/CHPclim_V1_NA.npz'
FileCorr_ratio = '/home/gut428/CHPclim_V1_NA_corratio.npz'
# FileClimo_tar = '/datastore/GLOBALWATER/CommonData/EMDNA_new/PBCOR_V1.0/WorldClim_V2_NA.npz'
# FileCorr_ratio = '/home/gut428/WorldClim_V2_NA_corratio.npz'
# FileClimo_tar = '/datastore/GLOBALWATER/CommonData/EMDNA_new/PBCOR_V1.0/CHELSA_V12_NA.npz'
# FileCorr_ratio = '/home/gut428/CHELSA_V12_NA_corratio.npz'

# meshed lat/lon of the target region
ncfid = nc.Dataset(FileGridInfo)
lattarm = ncfid.variables['latitude'][:].data
lattarm = np.flipud(lattarm)
lontarm = ncfid.variables['longitude'][:].data
ncfid.close()
lontar = lontarm[0, :]
lattar = lattarm[:, 0]
nrows = len(lattar)
ncols = len(lontar)

# grid infomation
ncfid = nc.Dataset(FileGridInfo)
gridlat = ncfid.variables['latitude'][:].data
gridlon = ncfid.variables['longitude'][:].data
gridele = ncfid.variables['elev'][:].data
gridgns = ncfid.variables['gradient_n_s'][:].data
gridgwe = ncfid.variables['gradient_w_e'][:].data
mask = ncfid.variables['mask'][:].data  # 1: grids to be considered; the other values: invalid grids
ncfid.close()

gridinfo = np.zeros([nrows, ncols, 6])
gridinfo[:, :, 0] = 1
gridinfo[:, :, 1] = gridlat
gridinfo[:, :, 2] = gridlon
gridinfo[:, :, 3] = gridele
gridinfo[:, :, 4] = gridgns
gridinfo[:, :, 5] = gridgwe
del gridlat, gridlon, gridele, gridgns, gridgwe

# calculate correction ratio for each month
if os.path.isfile(FileClimo_tar):
    datatemp = np.load(FileClimo_tar)
    pclim_tar=datatemp['pclim_tar']
    pclim_num=datatemp['pclim_num']
    pclim_tar[pclim_num < 2] = np.nan
    del datatemp
else:
    datatemp = nc.Dataset(FileClimo)
    pclim = datatemp['corr_P_monthly'][:].data
    latclim = datatemp['lat'][:].data
    lonclim = datatemp['lon'][:].data

    pclim_tar = np.nan * np.zeros([nrows, ncols, 12])
    pclim_num = np.nan * np.zeros([nrows, ncols])
    res = 0.1/2 # half resolution
    for rr in range(nrows):
        if np.mod(rr,50)==0:
            print(rr, nrows)
        for cc in range(ncols):
            if mask[rr, cc] != 1:
                continue
            latrc = gridinfo[rr, cc, 1]
            lonrc = gridinfo[rr, cc, 1]
            lonind = np.where((lonclim >= lonrc-res) & (lonclim <= lonrc+res))[0]
            latind = np.where((latclim >= latrc - res) & (latclim <= latrc + res))[0]
            pclim_rc = pclim[:, latind[0]:latind[1]+1, lonind[0]:lonind[1]+1]
            pclim_num[rr, cc] = np.sum(~np.isnan(pclim_rc[0, :, :]))
            for m in range(12):
                pclim_tar[rr, cc, m] = np.nanmean(pclim_rc[m, :, :])
    pclim_tar = np.flipud(pclim_tar)
    pclim_num = np.flipud(pclim_num)
    np.savez_compressed(FileClimo_tar,pclim_tar=pclim_tar,pclim_num=pclim_num,lattar=lattar,lontar=lontar)

if not os.path.isfile(FileCorr_ratio):
    corr_ratio = np.nan * np.zeros([nrows, ncols, 12])
    for month in range(12):
        print(month)
        data_oi = np.nan * np.zeros([nrows, ncols, 40])
        flag = 0
        for year in range(1979, 2019):
            # print(year)
            fileoi = path_oi + '/oimerge_prcp' + str(year * 100 + month + 1) + '.npz'
            datatemp = np.load(fileoi)
            data_oi[:,:,flag] = np.sum(datatemp['oi_value'],axis=2)
            flag = flag + 1
        data_oim = np.mean(data_oi, axis=2)
        ratiom = pclim_tar[:,:,month] / data_oim
        ratiom[data_oim == 0] = 0
        corr_ratio[:, :, month] = ratiom
    mask1 = np.flipud(mask)
    mask2 = corr_ratio[:, :, 0].copy()
    for rr in range(nrows):
        for cc in range(ncols):
            if mask1[rr,cc] == 1 and np.isnan(mask2[rr,cc]):
                rd = 1
                flag = 0
                while flag == 0:
                    drc = mask2[rr-rd:rr+rd+1, cc-rd:cc+rd+1]
                    if np.sum(~np.isnan(drc))>0:
                        flag = 1
                    else:
                        rd = rd + 1
                for m in range(12):
                    drcm = corr_ratio[rr-rd:rr+rd+1,cc-rd:cc+rd+1,m]
                    drcm=drcm[~np.isnan(mask2[rr-rd:rr+rd+1,cc-rd:cc+rd+1])]
                    corr_ratio[rr, cc, m] = np.nanmedian(drcm)
    np.savez_compressed(FileCorr_ratio,corr_ratio=corr_ratio,lattar=lattar,lontar=lontar)