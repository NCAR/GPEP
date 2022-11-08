# optimal interpolation merging
# merge background (reanalysis) and observation (regression estimates)
# computation time:
# oi-merge for stations: 3 variables adn 12 months = 36 jobs. ~0.5 hour per job
# oi-merge for grids: 3 variables adn 12 months = 36 jobs. ~0.5 hour per job

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
# month = 1

print(vars,month)

########################################################################################################################

# basic settings
weightmode = 'BMA' # method used to merge different reanalysis products
# vars = ['prcp', 'tmean', 'trange']

# ### Local Mac settings
# # input files/paths
# FileGridInfo = '/Users/localuser/Research/EMDNA/basicinfo/gridinfo_whole.nc'
# path_bac = '/Users/localuser/Research/EMDNA/merge' # data that will be used as background
# path_obs = '/Users/localuser/Research/EMDNA/regression' # data that will be used as observation
# near_file_GMET = '/Users/localuser/Research/EMDNA/regression/weight_nearstn.npz' # near station of stations/grids
# file_mask = './DEM/NA_DEM_010deg_trim.mat'
# FileStnInfo = '/Users/localuser/Research/EMDNA/basicinfo/stnlist_whole.txt'
# gmet_stndatafile = '/Users/localuser/Research/EMDNA/stndata_whole.npz'
#
# # output files/paths (can also be used as inputs once generated)
# path_oimerge = '/Users/localuser/Research/EMDNA/oimerge'
# ### Local Mac settings


### Plato settings
# input files/paths
FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/gridinfo_whole.nc'
FileStnInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_aftercheck.npz'
path_bac = '/datastore/GLOBALWATER/CommonData/EMDNA_new/ReanalysisCorrMerge/GWRLSBMA_merge'
path_obs = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck'
near_file_GMET = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck/nearstn_catalog.npz'
file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA_new/DEM/NA_DEM_010deg_trim.mat'

# output files/paths (can also be used as inputs once generated)
path_oimerge = '/home/gut428/OImerge_prcp_loo_20200424'
### Plato settings

file_regression_stn = path_obs + '/regression_stn.npz'
file_corrmerge_stn = [''] * len(vars)
for i in range(len(vars)):
    if vars[i] == 'pop':
        file_corrmerge_stn[i] = path_bac + '/merge_stn_' + vars[i] + '_GWR_' + weightmode + '.npz'
    else:
        file_corrmerge_stn[i] = path_bac + '/mergecorr_stn_' + vars[i] + '_GWRLS_' + weightmode + '.npz'

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
lontar = lontarm[0, :]
lattar = lattarm[:, 0]
lontarm[np.isnan(mask)] = np.nan
lattarm[np.isnan(mask)] = np.nan

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
# OI-merging at station scale

# this function use nearby stations to get OI-merged estimate at the target station
# to save time, this part is run month by month, and then monthly files are combined to one
# this part should be revised if re-running is needed.

for v in range(len(vars)):
    filemerge_stn = path_oimerge + '/OImerge_stn_GWRBMA_' + vars[v] + '.npz'
    if os.path.isfile(filemerge_stn):
        continue
    print('OI merge at stations:', vars[v])

    # load station original observations
    datatemp = np.load(gmet_stndatafile)
    if vars[v] == 'pop':
        observation_stn = datatemp['prcp_stn']
        observation_stn[observation_stn > 0] = 1
    else:
        observation_stn = datatemp[vars[v]+'_stn']

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
        near_loc = datatemp['near_stn_prcpLoc']
        near_weight = datatemp['near_stn_prcpWeight']
    else:
        near_loc = datatemp['near_stn_tempLoc']
        near_weight = datatemp['near_stn_tempWeight']
    del datatemp

    # start OI merging
    oimerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)
    for m in range(month-1, month):
        print('month', m + 1)
        indm = (date_mm == m + 1)
        nday = sum(indm)

        filemerge_stnm = path_oimerge + '/OImerge_stn_GWRBMA_' + vars[v] + '_' + str(m+1) + '.npz'
        if not os.path.isfile(filemerge_stnm):
            # use optimal interpolation to get OI-merged estimate at each station points
            for i in range(nstn):
                if np.mod(i,100)==0:
                    print('station',i,nstn)
                if np.isnan(observation_stn[i, 0]):
                    continue

                near_loci = near_loc[i, :]
                near_loci = near_loci[near_loci > -1]

                b_tar = reafinal_stn[i, indm]
                o_tar = regression_stn[i, indm]
                b_near = reafinal_stn[near_loci,:][:, indm]
                o_near = regression_stn[near_loci,:][:, indm]

                near_err_b = b_near - observation_stn[near_loci,:][:, indm]
                near_err_o = o_near - observation_stn[near_loci,:][:, indm]
                # tar_err_b = b_tar - observation_stn[i, indm]  # assume we know observation at the target point
                # in reality, we don't know target observation
                tar_err_b = np.zeros(nday)
                near_weighti = near_weight[i, 0:len(near_loci)]
                near_weighti = near_weighti / np.sum(near_weighti)
                for id in range(nday):
                    tar_err_b[id] = np.sum(near_err_b[:, id] * near_weighti)

                # delete possible nan values
                indnan = ~np.isnan(tar_err_b + np.sum(near_err_b,axis=0) + np.sum(near_err_o,axis=0) )
                weight = OImerge(tar_err_b[indnan], near_err_b[:, indnan], near_err_o[:, indnan], eye_o=0)
                if np.any(np.isnan(weight)) or np.any(abs(weight) > 2):
                    weight = near_weighti

                diff = o_near - b_near
                merge_est = b_tar.copy()
                for id in range(nday):
                    merge_est[id] = merge_est[id] + np.dot(weight, diff[:, id])
                oimerge_stn[i, indm] = merge_est

            if vars[v] == 'pop':
                oimerge_stn[oimerge_stn < 0] = 0
                oimerge_stn[oimerge_stn > 1] = 1
            if vars[v] == 'prcp':
                oimerge_stn[oimerge_stn < 0] = 0
            if vars[v] == 'trange':
                oimerge_stn = np.abs(oimerge_stn)
            np.savez_compressed(filemerge_stnm, oimerge_stn=oimerge_stn)


# combine all months
for v in range(len(vars)):
    flag = 1
    for m in range(12):
        filemerge_stnm = path_oimerge + '/OImerge_stn_GWRBMA_' + vars[v] + '_' + str(m+1) + '.npz'
        if not os.path.isfile(filemerge_stnm):
            flag = 0
    if flag==1:
        print('Combine station OI merging data for 12 months')
        oimerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)
        filemerge_stn = path_oimerge + '/OImerge_stn_GWRBMA_' + vars[v] + '.npz'
        if not os.path.isfile(filemerge_stn):
            for m in range(12):
                indm = (date_mm == m + 1)
                nday = sum(indm)
                filemerge_stnm = path_oimerge + '/OImerge_stn_GWRBMA_' + vars[v] + '_' + str(m + 1) + '.npz'
                datatemp = np.load(filemerge_stnm)
                oimerge_stnm = datatemp['oimerge_stn']
                oimerge_stn[:, indm] = oimerge_stnm[:, indm]
                del oimerge_stnm
            np.savez_compressed(filemerge_stn, oimerge_stn=oimerge_stn, stninfo=stninfo, date_ymd=date_ymd, stnID=stnID)

oup  # or output np.vstack((KGE, r, alpha, beta))