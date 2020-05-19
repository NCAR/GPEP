# merge background (reanalysis) and observation (regression estimates)

import numpy as np
import auxiliary as au
from matplotlib import pyplot as plt
from scipy import io
import os
import sys
import h5py
import time
import random
import datetime
from optimal_interpolation import OImerge
import calendar


def box_cox_transform(data, exp=0.25):
    return (data ** exp - 1) / exp


def box_cox_recover(data, exp=0.25):
    dataout = (data * exp + 1) ** (1 / exp)
    dataout[data < -1 / exp] = 0
    return dataout


# basic settings
# y1 = int(sys.argv[3])
# y2 = int(sys.argv[4])
# year = [y1, y2]
year = [1980, 1980]

weightmode = 'RMSE'
vars = ['prcp', 'tmean', 'trange']
# "Gaussian": prcp will be transformed into normal distributions; "Actual": actual space
prcp_space = 'Actual'
hwsize = 2  # 5X5 space window used to support estimation at the center grid
lontar = np.arange(-180 + 0.05, -50, 0.1)
lattar = np.arange(85 - 0.05, 5, -0.1)

# output path of merged data
path_merge = '/Users/localuser/GMET/merge'

# path of merged reanalysis
# path_bac = '/datastore/GLOBALWATER/CommonData/EMDNA/ReanalysisCorrMerge/Reanalysis_merge'
path_bac = '/Users/localuser/GMET/merge'

# path of gridded observations
# path_obs = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout'
path_obs = '/Users/localuser/GMET/merge'

# mask file
# file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA/DEM/NA_DEM_010deg_trim.mat'
file_mask = './DEM/NA_DEM_010deg_trim.mat'
mask = io.loadmat(file_mask)
mask = mask['DEM']
mask[~np.isnan(mask)] = 1  # 1: valid pixels
nrows, ncols = np.shape(mask)

# start
for y in range(year[0], year[1] + 1):
    for m in range(12):
        print('merging date:', y * 100 + m + 1)
        filemerge = path_merge + '/OImerge_' + str(y * 100 + m + 1) + '.npz'
        if os.path.isfile(filemerge):
            print('file exists')
            continue

        date_cal_start = y * 10000 + (m + 1) * 100 + 1
        date_cal_end = y * 10000 + (m + 1) * 100 + calendar.monthrange(y, m + 1)[1]
        datestr = str(date_cal_start) + '-' + str(date_cal_end)
        oimerge_data = dict()
        for var in vars:
            print('variable', var)
            file_bac = path_bac + '/mergedata_' + var + '_' + str(y * 100 + m + 1) + weightmode + '.npz'
            file_obs = path_obs + '/output_' + datestr + '.npz'
            file_obserr = path_obs + '/output_realerror_' + datestr + '.npz'

            # load background data (value and error)
            datatemp = np.load(file_bac)
            v_bac = datatemp['merge_data']
            if var == 'tmean' or var == 'trange':
                e_bac = datatemp['merge_error_raw']  # error in actual space
            elif var == 'prcp':
                if prcp_space == 'Actual':
                    e_bac = datatemp['merge_error_raw']
                elif prcp_space == 'Gaussian':
                    if weightmode == 'BMA':
                        e_bac = datatemp['merge_error_bc']
                        v_bac = box_cox_transform(v_bac)
                    elif weightmode == 'RMSE':
                        e_bac = datatemp['merge_error_raw']
                        truth = v_bac - e_bac
                        v_bac = box_cox_transform(v_bac)
                        e_bac = v_bac - box_cox_transform(truth)
            else:
                sys.exit('Unknown variable')
            del datatemp

            # load observation data (value and error)
            datatemp = np.load(file_obs)
            datatemp2 = np.load(file_obserr)
            if var == 'tmean' or var == 'trange':
                v_obs = datatemp[var]
                # e_obs = datatemp[var + '_err']
                e_obs = datatemp2[var + '_realerr']
            elif var == 'prcp':
                if prcp_space == 'Actual':
                    v_obs = datatemp['pcp_raw']
                    # e_obs = datatemp['pcp_err_raw']
                    e_obs = datatemp2['pcp_realerr_raw']
                elif prcp_space == 'Gaussian':
                    v_obs = datatemp['pcp_bc']
                    # e_obs = datatemp['pcp_err_bc']
                    e_obs = datatemp['pcp_realerr_bc']
            else:
                sys.exit('Unknown variable')
            del datatemp
            v_obs = np.flipud(v_obs)
            e_obs = np.flipud(e_obs)

            ntimes = np.shape(v_bac)[2]

            # produce OI merged estimateds
            oimerge_data[var] = np.nan * np.zeros(np.shape(v_bac), dtype=np.float32)
            weightsum = np.nan * np.zeros([nrows,ncols], dtype=np.float32)
            for r in range(nrows):
                if np.mod(r,50)==0:
                    print(r)
                for c in range(ncols):
                    if mask[r, c] != 1:
                        continue
                    if r < hwsize:
                        rse = [0, r + hwsize + 1]
                    elif r >= nrows - hwsize:
                        rse = [r - hwsize, nrows]
                    else:
                        rse = [r - hwsize, r + hwsize + 1]
                    if c < hwsize:
                        cse = [0, c + hwsize + 1]
                    elif c >= ncols - hwsize:
                        cse = [c - hwsize, ncols]
                    else:
                        cse = [c - hwsize, c + hwsize + 1]
                    snum = (rse[1] - rse[0]) * (cse[1] - cse[0])
                    tar_err_b = e_bac[r, c, :]
                    near_err_b = e_bac[rse[0]:rse[1], cse[0]:cse[1], :]
                    near_err_o = e_obs[rse[0]:rse[1], cse[0]:cse[1], :]
                    near_err_b = np.reshape(near_err_b, [snum, ntimes])
                    near_err_o = np.reshape(near_err_o, [snum, ntimes])
                    indnan = np.isnan(near_err_o[:, 0]) | np.isnan(near_err_b[:, 0])
                    near_err_b = near_err_b[~indnan, :]
                    near_err_o = near_err_o[~indnan, :]

                    weight = OImerge(tar_err_b, near_err_b, near_err_o)
                    # check the effectiveness of the weight since weight is not correct for very few cases
                    if np.isnan(weight[0]) or np.any(abs(weight)>2):
                        pass

                    diff = v_obs[rse[0]:rse[1], cse[0]:cse[1], :] - v_bac[rse[0]:rse[1], cse[0]:cse[1], :]
                    diff = np.reshape(diff, [snum, ntimes])
                    diff = diff[~indnan, :]
                    merge_est = v_bac[r, c, :]
                    for i in range(ntimes):
                        merge_est[i] = merge_est[i] + np.dot(weight, diff[:, i])

                    oimerge_data[var][r, c, :] = merge_est
                    weightsum[r,c] = np.sum(weight)
        if var == 'prcp' and prcp_space == 'Gaussian':
            oimerge_data[var] = box_cox_recover(oimerge_data[var])
    np.savez_compressed(filemerge, pcp=oimerge_data['prcp'], tmean=oimerge_data['tmean'], trange=oimerge_data['trange'])
