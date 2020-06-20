# optimal interpolation merging
# merge background (reanalysis) and observation (regression estimates)

import numpy as np
from scipy import io
import os
import netCDF4 as nc
from optimal_interpolation import OImerge
from auxiliary_merge import *

########################################################################################################################

# time periods and methods
vars = sys.argv[1]
vars = [vars]
month = int(sys.argv[2])

# vars = 'prcp'
# vars = [vars]
# month = 2

print(vars,month)

########################################################################################################################

# basic settings
weightmode = 'BMA' # method used to merge different reanalysis products
# vars = ['prcp', 'tmean', 'trange']

# ### Local Mac settings
# # input files/paths
# path_bac = '/Users/localuser/Research/EMDNA/merge' # data that will be used as background
# path_obs = '/Users/localuser/Research/EMDNA/regression' # data that will be used as observation
# near_file_GMET = '/Users/localuser/Research/EMDNA/regression/weight_nearstn.npz' # near station of stations/grids
# file_mask = './DEM/NA_DEM_010deg_trim.mat'
# FileStnInfo = '/Users/localuser/GMET/pyGMET_NA/stnlist_whole.txt'
# gmet_stndatafile = '/Users/localuser/Research/EMDNA/stndata_whole.npz'
#
# # output files/paths (can also be used as inputs once generated)
# path_oimerge = '/Users/localuser/Research/EMDNA/oimerge'
#
# ### Local Mac settings


### Plato settings
# input files/paths
FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/gridinfo_whole.nc'
FileStnInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_aftercheck.npz'
path_bac = '/datastore/GLOBALWATER/CommonData/EMDNA_new/ReanalysisCorrMerge/Reanalysis_merge'
path_obs = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck'
near_file_GMET = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck/nearstn_catalog.npz'
file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA_new/DEM/NA_DEM_010deg_trim.mat'

# output files/paths (can also be used as inputs once generated)
path_oimerge = '/home/gut428/OImerge'
### Plato settings

file_regression_stn = path_obs + '/regression_stn.npz'
file_corrmerge_stn = [''] * len(vars)
for i in range(len(vars)):
    if vars[i] == 'pop':
        file_corrmerge_stn[i] = path_bac + '/merge_stn_' + vars[i] + '_GWR_' + weightmode + '.npz'
    else:
        file_corrmerge_stn[i] = path_bac + '/mergecorr_stn_' + vars[i] + '_GWRQM_' + weightmode + '.npz'

########################################################################################################################

# basic processing
mask = io.loadmat(file_mask)
mask = mask['DEM']
mask[~np.isnan(mask)] = 1  # 1: valid pixels
nrows, ncols = np.shape(mask)

# meshed lat/lon of the target region
ncfid = nc.Dataset(FileGridInfo)
lattarm = ncfid.variables['latitude'][:].data
lattarm = np.flipud(lattarm)
lontarm = ncfid.variables['longitude'][:].data
ncfid.close()
lontarm[np.isnan(mask)] = np.nan
lattarm[np.isnan(mask)] = np.nan
lontar = lontarm[0, :]
lattar = lattarm[:, 0]

# load observations for all stations
datatemp = np.load(gmet_stndatafile)
stndata = datatemp['prcp_stn']
stninfo = datatemp['stninfo']
stnID = datatemp['stnID']
date_ymd = datatemp['date_ymd']
nstn, ntimes = np.shape(stndata)
del datatemp
date_yyyy = (date_ymd/10000).astype(int)
date_mm = (np.mod(date_ymd, 10000)/100).astype(int)

########################################################################################################################
# OI-merging at grid scale

for v in range(len(vars)):
    print('OI merge at grids:', vars[v])

    # load station original observations
    datatemp = np.load(gmet_stndatafile)
    if vars[v] == 'pop':
        observation_stn = datatemp['prcp_stn']
        observation_stn[observation_stn > 0] = 1
    else:
        observation_stn = datatemp[vars[v] + '_stn']
    del datatemp

    # load station regression estimates (obs)
    datatemp = np.load(file_regression_stn)
    regression_stn = datatemp[vars[v]]
    del datatemp

    # load corrected/merged reanalysis data at all station points (those are totally independent with station observations)
    datatemp = np.load(file_corrmerge_stn[v])
    reafinal_stn = datatemp['reamerge_stn']
    nstn, ntimes = np.shape(reafinal_stn)
    del datatemp

    # load near station information
    datatemp = np.load(near_file_GMET)
    if vars[v] == 'prcp' or vars[v] == 'pop':
        near_loc = datatemp['near_grid_prcpLoc']
        near_weight = datatemp['near_grid_prcpWeight']
        near_dist = datatemp['near_grid_prcpDist']
    else:
        near_loc = datatemp['near_grid_tempLoc']
        near_weight = datatemp['near_grid_tempWeight']
        near_dist = datatemp['near_grid_tempDist']
    near_loc = np.flipud(near_loc)
    near_weight = np.flipud(near_weight)
    near_dist = np.flipud(near_dist)
    del datatemp

    # start OI merging
    for m in range(month-1, month):
        print('month', m + 1)
        indm = (date_mm == m + 1)
        nday = sum(indm)
        datem = date_yyyy[indm]

        # load gridded merged reanalysis data for all years
        print('load gridded merged data')
        reagrid_value = np.nan * np.zeros([nrows, ncols, nday], dtype=np.float32)
        reagrid_error = np.nan * np.zeros([nrows, ncols, nday], dtype=np.float32)
        for y in range(1979, 2019):
            indym = datem == y
            filey = path_bac + '/bmamerge_' + vars[v] + '_' + str(y*100+m+1) + '.npz'
            datatemp = np.load(filey)
            reagrid_value[:, :, indym] = datatemp['bma_data']
            reagrid_error[:, :, indym] = datatemp['bma_error']
            del datatemp

        # calculate OI-merging weights for every grids
        print('calculate OI merging weights')
        file_oiweight = path_oimerge + '/oiweight_' + vars[v] + '_month_' + str(m+1) + '.npz'
        if os.path.isfile(file_oiweight):
            datatemp = np.load(file_oiweight)
            oiweight = datatemp['oiweight']
            del datatemp
        else:
            oiweight = np.nan * np.zeros([nrows, ncols, np.shape(near_loc)[2]], dtype=np.float32)
            for r in range(nrows):
                if np.mod(r,10)==0:
                    print(r)
                for c in range(ncols):
                    if np.isnan(mask[r, c]):
                        continue
                    near_loci = near_loc[r, c, :]
                    near_loci = near_loci[near_loci > -1]

                    b_near = reafinal_stn[near_loci, :][:, indm]
                    o_near = regression_stn[near_loci, :][:, indm]
                    # this error is from weighted mean. if using nearest neighbor to obtain gridded error, this weight will be more similar to stn-OI
                    tar_err_b = reagrid_error[r, c, :]
                    near_err_b = b_near - observation_stn[near_loci, :][:, indm]
                    near_err_o = o_near - observation_stn[near_loci, :][:, indm]

                    # delete possible nan values
                    induse = ~np.isnan(tar_err_b + np.sum(near_err_b, axis=0) + np.sum(near_err_o, axis=0))
                    weight = OImerge(tar_err_b[induse], near_err_b[:, induse], near_err_o[:, induse], eye_o=0)

                    oiweight[r, c, 0:len(weight)] = weight
            np.savez_compressed(file_oiweight, oiweight=oiweight, lattar=lattar, lontar=lontar)