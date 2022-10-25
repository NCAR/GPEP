# perform locally weighted regression for all stations using leave-one-out
# regression estimates can support further screening of stations and Optimal Interpolation
# calculation time: 11 hours when submitted as 40 Jobs (40 years)


import numpy as np
import regression as reg
import os
import sys

########################################################################################################################

# setting: parameters for transforming temp to approximate normal distribution
trans_mode = 'none'  # box-cox or power-law or none
trans_exp_daily = np.nan

# infiles
gmet_stndatafile = './Andy_test_case/stndata_whole.npz'
file_nearstn = './Andy_test_case/nearstn_catalog.npz'

# outfiles
outpath = './Andy_test_case/raw_regression'
os.makedirs(outpath, exist_ok=True)

filereg_all = outpath + '/regression_stn.npz'  # regression error at station points

########################################################################################################################

# load station data
if not os.path.isfile(gmet_stndatafile):
    print(gmet_stndatafile,'does not exist')
    sys.exit()
else:
    print('load station data')
    datatemp=np.load(gmet_stndatafile)
    stninfo=datatemp['stninfo']
    stnID=datatemp['stnID']
    nstn = len(stnID)

    prcp_stn_daily = datatemp['prcp_stn'][:]
    tmean_stn_daily = datatemp['tmean_stn'][:]
    trange_stn_daily = datatemp['trange_stn'][:]

########################################################################################################################

# find neighboring stations and calculate distance-based weights
if not os.path.isfile(file_nearstn):
    print(file_nearstn, 'does not exist')
    sys.exit()
else:
    print('file_nearstn exists. loading ...')
    with np.load(file_nearstn) as datatemp:
        near_stn_prcpLoc = datatemp['near_stn_prcpLoc']
        near_stn_prcpWeight = datatemp['near_stn_prcpWeight']
        near_stn_tempLoc = datatemp['near_stn_tempLoc']
        near_stn_tempWeight = datatemp['near_stn_tempWeight']
    del datatemp

########################################################################################################################

if not os.path.isfile(filereg_all):
    print('Estimate daily regression error at station points')
    prcp_err_stn_daily, tmean_err_stn_daily, trange_err_stn_daily, pop_err_stn_daily = \
        reg.station_error(prcp_stn_daily, tmean_stn_daily, trange_stn_daily, stninfo, near_stn_prcpLoc,
                          near_stn_prcpWeight, near_stn_tempLoc, near_stn_tempWeight, trans_exp_daily,
                          trans_mode, 10)

    prcp_reg_stn = prcp_err_stn_daily + prcp_stn_daily
    prcp_reg_stn[prcp_reg_stn < 0] = 0
    tmean_reg_stn = tmean_err_stn_daily + tmean_stn_daily
    trange_reg_stn = trange_err_stn_daily + trange_stn_daily

    prcp_stn_daily[prcp_stn_daily>0] = 1
    pop_reg_stn = pop_err_stn_daily + prcp_stn_daily

    np.savez_compressed(filereg_all, prcp=prcp_reg_stn, tmean=tmean_reg_stn,
                        trange=trange_reg_stn, pop=pop_reg_stn)