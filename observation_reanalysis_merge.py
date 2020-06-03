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
y1 = int(sys.argv[1])
y2 = int(sys.argv[2])
vars = sys.argv[3]
vars = [vars]
month = int(sys.argv[4])
year = [y1, y2]
# year = [1980, 1980]
print('year:',year)

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
path_bac = '/home/gut428/ReanalysisCorrMerge/Reanalysis_merge'
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

# get the error of the merged data for all stations and extrapolate to all grids
# this function use nearby stations to get OI-merged estimate at the target station
for v in range(len(vars)):
    print('OI merge at stations:', vars[v])
    # filemerge_stn = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '.npz'
    # if os.path.isfile(filemerge_stn):
    #     continue

    # load station original observations
    datatemp = np.load(gmet_stndatafile)
    observation_stn = datatemp[vars[v]+'_stn']

    # load station regression estimates (obs)
    datatemp = np.load(file_regression_stn)
    regression_stn = datatemp[vars[v]]
    del datatemp

    # load corrected/merged reanalysis data at all station points (those are totally independent with station observations)
    # and find the best choice
    datatemp = np.load(file_corrmerge_stn[v])
    reamerge_stn_all = datatemp['reamerge_stn']
    reacorr_stn_all = datatemp['reacorr_stn']
    reanum, nstn, ntimes = np.shape(reacorr_stn_all)
    del datatemp
    reafinal_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)
    for m in range(12):
        indm =  (date_number['mm'] == m + 1)
        rearmse = np.zeros([nstn, reanum + 1])
        rearmse[:, 0] = calmetric(reamerge_stn_all[:, indm], observation_stn[:, indm], metname='RMSE')
        for i in range(reanum):
            rearmse[:, i + 1] = calmetric(reacorr_stn_all[i, :, indm].T, observation_stn[:, indm], metname='RMSE')
        bestchoice = np.argmin(rearmse, axis=1)
        for i in range(nstn):
            if bestchoice[i] > 0:
                reafinal_stn[i, indm] = reacorr_stn_all[bestchoice[i] - 1, i, indm]
            else:
                reafinal_stn[i, indm] = reamerge_stn_all[i, indm]

    # load near station information
    datatemp = np.load(near_file_GMET)
    if vars[v] == 'prcp':
        near_loc = datatemp['near_stn_prcpLoc']
        near_weight = datatemp['near_stn_prcpWeight']
    else:
        near_loc = datatemp['near_stn_tempLoc']
        near_weight = datatemp['near_stn_tempWeight']
    del datatemp

    # start OI merging
    oimerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.nan)

    for m in range(month, month+1):
        print('month', m + 1)
        indm = (date_number['mm'] == m + 1)
        nday = sum(indm)

        filemerge_stnm = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '_' + str(m+1) + '.npz'
        if not os.path.isfile(filemerge_stnm):
            # use optimal interpolation to get OI-merged estimate at each station points
            for i in range(nstn):
                if np.mod(i,2000)==0:
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

            np.savez_compressed(filemerge_stnm, oimerge_stn=oimerge_stn, reafinal_stn=reafinal_stn)


# combine all months
for v in range(len(vars)):
    flag = 1
    for m in range(12):
        print('month', m + 1)
        filemerge_stnm = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '_' + str(m+1) + '.npz'
        if not os.path.isfile(filemerge_stnm):
            flag = 0

    if flag==1:
        oimerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.nan)
        filemerge_stn = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '.npz'
        for m in range(12):
            indm = (date_number['mm'] == m + 1)
            nday = sum(indm)
            filemerge_stnm = path_oimerge + '/OImerge_stn_GWRQMBMA_' + vars[v] + '_' + str(m + 1) + '.npz'
            datatemp = np.load(filemerge_stnm)
            oimerge_stnm = datatemp['oimerge_stn']
            oimerge_stn[:, indm] = oimerge_stnm[:, indm]
        np.savez_compressed(filemerge_stn, oimerge_stn=oimerge_stn, reafinal_stn=reafinal_stn)

########################################################################################################################

# get the error of the merged data for all stations and extrapolate to all grids
# this function use nearby grids to get OI-merged estimate at the target station (this is closer to actual gridded OI)

# datatemp = np.load(file_corrmerge_stn[0])
# reacorr_stn_all = datatemp['reacorr_stn']
# reanum, nstn, ntimes = np.shape(reacorr_stn_all)
# del datatemp, reacorr_stn_all
#
# # find nearby grids for each station
# nearrowcol = np.zeros([nstn, 2], dtype=int)
# for i in range(nstn):
#     rowi = np.argmin(abs(lattar-stninfo[i, 1]))
#     coli = np.argmin(abs(lontar-stninfo[i, 2]))
#     nearrowcol[i, 0] = rowi
#     nearrowcol[i, 1] = coli
#
# # start
# for v in range(len(vars)):
#     print('OI merge at stations:', vars[v])
#     filemerge_stn = path_oimerge + '/OImerge_stn_' + vars[v] + '_neargrid.npz'
#     if os.path.isfile(filemerge_stn):
#         continue
#
#     # load corrected/merged reanalysis data at all station points (those are totally independent with station observations)
#     datatemp = np.load(file_corrmerge_stn[v])
#     reamerge_stn_all = datatemp['reamerge_stn']
#     reacorr_stn_all = datatemp['reacorr_stn']
#     reanum, nstn, ntimes = np.shape(reacorr_stn_all)
#     del datatemp
#
#     oimerge_stn = np.zeros([nstn, ntimes])
#
#     for y in range(1979, 2019):
#         for m in range(12):
#             print('merging error: date:', y * 100 + m + 1)
#             indym = (date_number['yyyy'] == y) & (date_number['mm'] == m + 1)
#             nday = sum(indym)
#             oimerge_stnym = np.zeros([nstn, nday])
#
#             # load station regression estimates at all stations (locally weighted regression)
#             date_cal_start = y * 10000 + (m + 1) * 100 + 1
#             date_cal_end = y * 10000 + (m + 1) * 100 + calendar.monthrange(y, m + 1)[1]
#             datestr = str(date_cal_start) + '-' + str(date_cal_end)
#             file1 = path_obs + '/stndata_' + datestr + '.npz'
#             datatemp1 = np.load(file1)
#             observation_stn = datatemp1[vars[v] + '_stn_daily']
#
#
#
#             file_bac = path_bac + '/mergedata_' + vars[v] + '_' + str(y * 100 + m + 1) + weightmode + '.npz'
#             file_obs = path_obs + '/output_' + datestr + '.npz'
#             file_obserr = path_obs + '/output_realerror_' + datestr + '.npz'
#
#
#             # get the best reanalysis among three original reanalysis datasets and one merged dataset for each station
#             reamerge_stn = reamerge_stn_all[:, indym]
#             reacorr_stn = reacorr_stn_all[:, :, indym]
#             rearmse = np.zeros([nstn, reanum+1])
#             rearmse[:, 0] = calmetric(reamerge_stn, observation_stn, metname='RMSE')
#             for i in range(3):
#                 rearmse[:, i + 1] = calmetric(reacorr_stn[i, :, :], observation_stn, metname='RMSE')
#             bestchoice = np.argmin(rearmse, axis=1)
#             reafinal_stn = reamerge_stn.copy()
#             for i in range(nstn):
#                 if bestchoice[i] > 0:
#                     reafinal_stn[i, :] = reacorr_stn[bestchoice[i] - 1, i, :]
#
#             # load background data (value and error)
#             datatemp = np.load(file_bac)
#             v_bac = datatemp['merge_data']
#             if vars[v] == 'tmean' or vars[v] == 'trange':
#                 e_bac = datatemp['merge_error_raw']  # error in actual space
#             elif vars[v] == 'prcp':
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
#
#             # load observation data (value and error)
#             datatemp = np.load(file_obs)
#             datatemp2 = np.load(file_obserr)
#             if vars[v] == 'tmean' or vars[v] == 'trange':
#                 v_obs = datatemp[vars[v]]
#                 # e_obs = datatemp[var + '_err']
#                 e_obs = datatemp2[vars[v] + '_realerr']
#             elif vars[v] == 'prcp':
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
#             # use optimal interpolation to get OI-merged estimate at each station points
#             for i in range(nstn):
#                 if not np.isnan(observation_stn[i, 0]):
#                     b_tar = reafinal_stn[i, :]
#                     tar_err_b = b_tar - observation_stn[i, :]
#
#                     b_near = v_bac[nearrowcol[i, 0]-hwsize:nearrowcol[i, 0]+hwsize+1,
#                              nearrowcol[i, 1]-hwsize:nearrowcol[i, 1]+hwsize+1, :]
#                     o_near = v_obs[nearrowcol[i, 0]-hwsize:nearrowcol[i, 0]+hwsize+1,
#                              nearrowcol[i, 1]-hwsize:nearrowcol[i, 1]+hwsize+1, :]
#                     near_err_b = e_bac[nearrowcol[i, 0]-hwsize:nearrowcol[i, 0]+hwsize+1,
#                              nearrowcol[i, 1]-hwsize:nearrowcol[i, 1]+hwsize+1, :]
#                     near_err_o = e_obs[nearrowcol[i, 0]-hwsize:nearrowcol[i, 0]+hwsize+1,
#                              nearrowcol[i, 1]-hwsize:nearrowcol[i, 1]+hwsize+1, :]
#                     b_near[hwsize, hwsize, :] = np.nan
#                     o_near[hwsize, hwsize, :] = np.nan
#                     near_err_b[hwsize, hwsize, :] = np.nan
#                     near_err_o[hwsize, hwsize, :] = np.nan
#                     snum=(hwsize*2+1)**2
#                     b_near=np.reshape(b_near,[snum, nday])
#                     o_near = np.reshape(o_near, [snum, nday])
#                     near_err_b = np.reshape(near_err_b, [snum, nday])
#                     near_err_o = np.reshape(near_err_o, [snum, nday])
#                     indnan = ~np.isnan(b_near[:,0])
#                     b_near=b_near[indnan,:]
#                     o_near = o_near[indnan, :]
#                     near_err_b = near_err_b[indnan, :]
#                     near_err_o = near_err_o[indnan, :]
#
#                     weight = OImerge(tar_err_b, near_err_b, near_err_o, 0)
#                     if np.any(np.isnan(weight)) or np.any(abs(weight) > 2):
#                         weight = np.ones(len(weight)) / len(weight)
#
#                     diff = o_near - b_near
#                     merge_est = b_tar.copy()
#                     for id in range(nday):
#                         merge_est[id] = merge_est[id] + np.dot(weight, diff[:, id])
#
#                     oimerge_stnym[i, :] = merge_est
#
#             oimerge_stn[:, indym] = oimerge_stnym
#
#     np.savez_compressed(filemerge_stn, oimerge_stn=oimerge_stn)


########################################################################################################################

# use optimal interpolation to obtain merged data for all grids
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
