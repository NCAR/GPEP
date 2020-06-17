
import numpy as np
from matplotlib import pyplot as plt
from scipy import io
import os
import sys
from optimal_interpolation import OImerge
import calendar
import auxiliary as au
from auxiliary_merge import *

########################################################################################################################

# time periods and methods
vars = sys.argv[1]
vars = [vars]
month = int(sys.argv[2])

# vars = 'prcp'
# vars = [vars]
# month = 2

########################################################################################################################

# basic settings
weightmode = 'BMA' # method used to merge different reanalysis products
# vars = ['prcp', 'tmean', 'trange']
hwsize = 2  # 5X5 space window used to support estimation at the center grid
lontar = np.arange(-180 + 0.05, -50, 0.1)
lattar = np.arange(85 - 0.05, 5, -0.1)

# "Gaussian": prcp will be transformed into normal distributions; "Actual": actual space
# "Gaussian" is not a good choice because station prcp regression using box-cox has large underestimation
prcp_space = 'Actual'

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
path_bac = '/datastore/GLOBALWATER/CommonData/EMDNA/ReanalysisCorrMerge/Reanalysis_merge'
path_obs = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout'
near_file_GMET = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout/weight.npz'
file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA/DEM/NA_DEM_010deg_trim.mat'
FileStnInfo = '/home/gut428/GMET/eCAI_EMDNA/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA/stndata_whole.npz'

# output files/paths (can also be used as inputs once generated)
path_oimerge = '/home/gut428/OImerge'
### Plato settings

file_regression_stn = path_obs + '/daily_regression_stn.npz'
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

# date
date_list, date_number = m_DateList(1979, 2018, 'ByYear')

# stninfo
stnID, stninfo = au.readstnlist(FileStnInfo)
nstn = len(stnID)

########################################################################################################################
# OI-merging at station scale

# this function use nearby stations to get OI-merged estimate at the target station
# to save time, this part is run month by month, and then monthly files are combined to one
# this part should be revised if re-running is needed.

for v in range(len(vars)):
    print('OI merge at stations:', vars[v])
    filemerge_stn = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '.npz'
    if os.path.isfile(filemerge_stn):
        continue

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
        indm = (date_number['mm'] == m + 1)
        nday = sum(indm)

        filemerge_stnm = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '_' + str(m+1) + '.npz'
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

                tar_err_b = b_tar - observation_stn[i, indm]
                near_err_b = b_near - observation_stn[near_loci,:][:, indm]
                near_err_o = o_near - observation_stn[near_loci,:][:, indm]

                # delete possible nan values
                indnan = ~np.isnan(tar_err_b + np.sum(near_err_b,axis=0) + np.sum(near_err_o,axis=0) )
                weight = OImerge(tar_err_b[indnan], near_err_b[:, indnan], near_err_o[:, indnan], eye_o=0)
                if np.any(np.isnan(weight)) or np.any(abs(weight) > 2):
                    weight = near_weight[i, 0:len(near_loci)]
                    weight = weight / np.sum(weight)

                diff = o_near - b_near
                merge_est = b_tar.copy()
                for id in range(nday):
                    merge_est[id] = merge_est[id] + np.dot(weight, diff[:, id])

                oimerge_stn[i, indm] = merge_est

            if vars[v] == 'pop':
                oimerge_stn[oimerge_stn < 0] = 0
                oimerge_stn[oimerge_stn > 1] = 1
            np.savez_compressed(filemerge_stnm, oimerge_stn=oimerge_stn)


# combine all months
for v in range(len(vars)):
    flag = 1
    for m in range(12):
        print('month', m + 1)
        filemerge_stnm = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '_' + str(m+1) + '.npz'
        if not os.path.isfile(filemerge_stnm):
            flag = 0
    if flag==1:
        oimerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)
        filemerge_stn = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '.npz'
        for m in range(12):
            indm = (date_number['mm'] == m + 1)
            nday = sum(indm)
            filemerge_stnm = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '_' + str(m + 1) + '.npz'
            datatemp = np.load(filemerge_stnm)
            oimerge_stnm = datatemp['oimerge_stn']
            oimerge_stn[:, indm] = oimerge_stnm[:, indm]
            del oimerge_stnm
        np.savez_compressed(filemerge_stn, oimerge_stn=oimerge_stn)