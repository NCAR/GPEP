# merge background (reanalysis) and observation (regression estimates)

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

# for v in range(len(vars)):
#     print('OI merge at stations:', vars[v])
#     # filemerge_stn = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '.npz'
#     # if os.path.isfile(filemerge_stn):
#     #     continue
#
#     # load station original observations
#     datatemp = np.load(gmet_stndatafile)
#     observation_stn = datatemp[vars[v]+'_stn']
#
#     # load station regression estimates (obs)
#     datatemp = np.load(file_regression_stn)
#     regression_stn = datatemp[vars[v]]
#     del datatemp
#
#     # load corrected/merged reanalysis data at all station points (those are totally independent with station observations)
#     datatemp = np.load(file_corrmerge_stn[v])
#     reafinal_stn = datatemp['reamerge_stn']
#     nstn, ntimes = np.shape(reafinal_stn)
#     del datatemp
#
#     # load near station information
#     datatemp = np.load(near_file_GMET)
#     if vars[v] == 'prcp':
#         near_loc = datatemp['near_stn_prcpLoc']
#         near_weight = datatemp['near_stn_prcpWeight']
#     else:
#         near_loc = datatemp['near_stn_tempLoc']
#         near_weight = datatemp['near_stn_tempWeight']
#     del datatemp
#
#     # start OI merging
#     oimerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)
#     for m in range(month-1, month):
#         print('month', m + 1)
#         indm = (date_number['mm'] == m + 1)
#         nday = sum(indm)
#
#         filemerge_stnm = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '_' + str(m+1) + '.npz'
#         if not os.path.isfile(filemerge_stnm):
#             # use optimal interpolation to get OI-merged estimate at each station points
#             for i in range(nstn):
#                 if np.mod(i,100)==0:
#                     print('station',i,nstn)
#
#                 if np.isnan(observation_stn[i, 0]):
#                     continue
#
#                 near_loci = near_loc[i, :]
#                 near_loci = near_loci[near_loci > -1]
#
#                 b_tar = reafinal_stn[i, indm]
#                 o_tar = regression_stn[i, indm]
#                 b_near = reafinal_stn[near_loci,:][:, indm]
#                 o_near = regression_stn[near_loci,:][:, indm]
#
#                 tar_err_b = b_tar - observation_stn[i, indm]
#                 near_err_b = b_near - observation_stn[near_loci,:][:, indm]
#                 near_err_o = o_near - observation_stn[near_loci,:][:, indm]
#
#                 # delete possible nan values
#                 indnan = ~np.isnan(tar_err_b + np.sum(near_err_b,axis=0) + np.sum(near_err_o,axis=0) )
#                 weight = OImerge(tar_err_b[indnan], near_err_b[:, indnan], near_err_o[:, indnan], eye_o=0)
#                 if np.any(np.isnan(weight)) or np.any(abs(weight) > 2):
#                     weight = near_weight[i, 0:len(near_loci)]
#                     weight = weight / np.sum(weight)
#
#                 diff = o_near - b_near
#                 merge_est = b_tar.copy()
#                 for id in range(nday):
#                     merge_est[id] = merge_est[id] + np.dot(weight, diff[:, id])
#
#                 oimerge_stn[i, indm] = merge_est
#
#             np.savez_compressed(filemerge_stnm, oimerge_stn=oimerge_stn)
#
#
# # combine all months
# for v in range(len(vars)):
#     flag = 1
#     for m in range(12):
#         print('month', m + 1)
#         filemerge_stnm = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '_' + str(m+1) + '.npz'
#         if not os.path.isfile(filemerge_stnm):
#             flag = 0
#     if flag==1:
#         oimerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)
#         filemerge_stn = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '.npz'
#         for m in range(12):
#             indm = (date_number['mm'] == m + 1)
#             nday = sum(indm)
#             filemerge_stnm = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '_' + str(m + 1) + '.npz'
#             datatemp = np.load(filemerge_stnm)
#             oimerge_stnm = datatemp['oimerge_stn']
#             oimerge_stn[:, indm] = oimerge_stnm[:, indm]
#         np.savez_compressed(filemerge_stn, oimerge_stn=oimerge_stn)

########################################################################################################################
# OI-merging at grid scale

for v in range(len(vars)):
    print('OI merge at stations:', vars[v])

    # load station original observations
    datatemp = np.load(gmet_stndatafile)
    observation_stn = datatemp[vars[v]+'_stn']
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
    if vars[v] == 'prcp':
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

    # load OI merged data at station points
    filemerge_stn = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '.npz'
    datatemp = np.load(filemerge_stn)
    oimerge_stn = datatemp['oimerge_stn']
    del datatemp

    # start OI merging
    for m in range(month-1, month):
        print('month', m + 1)
        indm = (date_number['mm'] == m + 1)
        nday = sum(indm)
        datem = date_number['yyyy'][indm]

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
                if np.mod(r,50)==0:
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
            np.savez_compressed(file_oiweight, oiweight=oiweight)

        # perform OI merging for all years
        print('perform OI merging')
        for y in range(1979, 2019):
            fileoi_ym = path_oimerge + '/oimerge_' + vars[v] + str(y*100+m+1) + '.npz'
            indym1 = datem == y
            ndayy = np.sum(indym1)
            indym2 = (date_number['mm'] == m + 1) & (date_number['yyyy'] == y)
            if os.path.isfile(fileoi_ym):
                continue

            # calculate OI value
            oi_value = np.nan * np.zeros([nrows, ncols, ndayy], dtype=np.float32)
            for r in range(nrows):
                for c in range(ncols):
                    if np.isnan(mask[r, c]):
                        continue
                    near_loci = near_loc[r, c, :]
                    near_loci = near_loci[near_loci > -1]

                    weight = oiweight[r, c, :]
                    weight = weight[~np.isnan(weight)]
                    if (np.any(np.isnan(weight))) or (np.any(abs(weight) > 2)): # exclude too large values
                        weight = near_weight[r, c, 0:len(weight)]
                        weight = weight / np.sum(weight)

                    b_tar = reagrid_value[r, c, indym1]
                    b_near = reafinal_stn[near_loci, :][:, indym2]
                    o_near = regression_stn[near_loci, :][:, indym2]

                    diff = o_near - b_near
                    merge_est = b_tar.copy()
                    for id in range(ndayy):
                        merge_est[id] = merge_est[id] + np.dot(weight, diff[:, id])
                    oi_value[r, c, :] = merge_est

            # calculate OI error (mean square error from nearby stations)
            oi_error_stn = (box_cox_transform(oimerge_stn[:, indym2]) - box_cox_transform(observation_stn[:, indym2])) ** 2
            oi_error_bc = extrapolation(oi_error_stn, near_loc, near_dist)

            oi_value = box_cox_transform(oi_value)
            np.savez_compressed(fileoi_ym, oi_value=oi_value, oi_error_bc=oi_error_bc)


# for y in range(year[0], year[1] + 1):
#     for m in range(12):
#         print('merging date:', y * 100 + m + 1)
#         filemerge = path_oimerge + '/OImerge_' + str(y * 100 + m + 1) + '.npz'
#         if os.path.isfile(filemerge):
#             print('file exists')
#             continue
#
#         date_cal_start = y * 10000 + (m + 1) * 100 + 1
#         date_cal_end = y * 10000 + (m + 1) * 100 + calendar.monthrange(y, m + 1)[1]
#         datestr = str(date_cal_start) + '-' + str(date_cal_end)
#         oimerge_data = dict()
#         for var in vars:
#             print('variable', var)
#             # due to the temporal limitation in data processing:
#             # the error in file_bac is mean error, and the error in file_bacerr is mean squared error
#             # the error in file_obs is mean squared error, and the error in file_obserr is mean error
#             file_bac = path_bac + '/mergedata_' + var + '_' + str(y * 100 + m + 1) + weightmode + '.npz'
#             file_bacerr = path_bac + '/mserror_' + var + '_' + str(y * 100 + m + 1) + weightmode + '.npz'
#             file_obs = path_obs + '/output_' + datestr + '.npz'
#             file_obserr = path_obs + '/output_realerror_' + datestr + '.npz'
#
#             # load background data (value and error)
#             datatemp = np.load(file_bac)
#             v_bac = datatemp['merge_data']
#             if var == 'tmean' or var == 'trange':
#                 e_bac = datatemp['merge_error_raw']  # error in actual space
#             elif var == 'prcp':
#                 if prcp_space == 'Actual':
#                     e_bac = datatemp['merge_error_raw']
#                 elif prcp_space == 'Gaussian':
#                     if weightmode == 'BMA':
#                         e_bac = datatemp['merge_error_bc']
#                         v_bac = box_cox_transform(v_bac)
#                     elif weightmode == 'RMSE':
#                         e_bac = datatemp['merge_error_raw']
#                         truth = v_bac - e_bac
#                         v_bac = box_cox_transform(v_bac)
#                         e_bac = v_bac - box_cox_transform(truth)
#             else:
#                 sys.exit('Unknown variable')
#             del datatemp
#
#             # load observation data (value and error)
#             datatemp = np.load(file_obs)
#             datatemp2 = np.load(file_obserr)
#             if var == 'tmean' or var == 'trange':
#                 v_obs = datatemp[var]
#                 # e_obs = datatemp[var + '_err']
#                 e_obs = datatemp2[var + '_realerr']
#             elif var == 'prcp':
#                 if prcp_space == 'Actual':
#                     v_obs = datatemp['pcp_raw']
#                     # e_obs = datatemp['pcp_err_raw']
#                     e_obs = datatemp2['pcp_realerr_raw']
#                 elif prcp_space == 'Gaussian':
#                     v_obs = datatemp['pcp_bc']
#                     # e_obs = datatemp['pcp_err_bc']
#                     e_obs = datatemp2['pcp_realerr_bc']
#             else:
#                 sys.exit('Unknown variable')
#             del datatemp
#             v_obs = np.flipud(v_obs)
#             e_obs = np.flipud(e_obs)
#
#             ntimes = np.shape(v_bac)[2]
#
#             # produce OI merged estimateds
#             oimerge_data[var] = np.nan * np.zeros(np.shape(v_bac), dtype=np.float32)
#             # weightsum = np.nan * np.zeros([nrows,ncols], dtype=np.float32)
#             for r in range(nrows):
#                 if np.mod(r, 50) == 0:
#                     print(r)
#                 for c in range(ncols):
#                     if mask[r, c] != 1:
#                         continue
#                     if r < hwsize:
#                         rse = [0, r + hwsize + 1]
#                     elif r >= nrows - hwsize:
#                         rse = [r - hwsize, nrows]
#                     else:
#                         rse = [r - hwsize, r + hwsize + 1]
#                     if c < hwsize:
#                         cse = [0, c + hwsize + 1]
#                     elif c >= ncols - hwsize:
#                         cse = [c - hwsize, ncols]
#                     else:
#                         cse = [c - hwsize, c + hwsize + 1]
#                     snum = (rse[1] - rse[0]) * (cse[1] - cse[0])
#                     tar_err_b = e_bac[r, c, :]
#                     near_err_b = e_bac[rse[0]:rse[1], cse[0]:cse[1], :]
#                     near_err_o = e_obs[rse[0]:rse[1], cse[0]:cse[1], :]
#                     near_err_b[r - rse[0], c - cse[0], :] = np.nan # tar grids should be excluded from nearby grids
#                     near_err_o[r - rse[0], c - cse[0], :] = np.nan
#                     near_err_b = np.reshape(near_err_b, [snum, ntimes])
#                     near_err_o = np.reshape(near_err_o, [snum, ntimes])
#                     indnan = np.isnan(near_err_o[:, 0]) | np.isnan(near_err_b[:, 0])
#                     near_err_b = near_err_b[~indnan, :]
#                     near_err_o = near_err_o[~indnan, :]
#
#                     weight = OImerge(tar_err_b, near_err_b, near_err_o, 0)
#                     # check the effectiveness of the weight since weight is not correct for very few cases
#                     # this is the disadvantage of using actual error estimation since there is large randomness
#                     if np.any(np.isnan(weight)) or np.any(abs(weight) > 2):
#                         weight0 = np.ones([rse[1] - rse[0], cse[1] - cse[0]])
#                         weight0[np.arange(rse[0], rse[1]) == r, np.arange(cse[0], cse[1]) == c] = 4
#                         weight = np.reshape(weight, snum)
#                         weight = weight[indnan]
#                         weight = weight / np.sum(weight)
#
#                     diff = v_obs[rse[0]:rse[1], cse[0]:cse[1], :] - v_bac[rse[0]:rse[1], cse[0]:cse[1], :]
#                     diff = np.reshape(diff, [snum, ntimes])
#                     diff = diff[~indnan, :]
#                     merge_est = v_bac[r, c, :].copy()
#                     for i in range(ntimes):
#                         merge_est[i] = merge_est[i] + np.dot(weight, diff[:, i])
#
#                     oimerge_data[var][r, c, :] = merge_est
#                     # weightsum[r, c] = np.sum(weight)
#             if var == 'prcp' and prcp_space == 'Gaussian':
#                 oimerge_data[var] = box_cox_recover(oimerge_data[var])
#     np.savez_compressed(filemerge, pcp=oimerge_data['prcp'], tmean=oimerge_data['tmean'], trange=oimerge_data['trange'],
#                         latitude=lattar, longitude=lontar)
