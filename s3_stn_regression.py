# perform locally weighted regression for all stations using leave-one-out
# regression estimates can support further screening of stations and Optimal Interpolation

import numpy as np
import netCDF4 as nc
import auxiliary as au
import regression as reg
import datetime as dt
import os
import sys
from auxiliary_merge import m_DateList


########################################################################################################################

year = int(sys.argv[1])
yearall = [1979, 2018] # this is used for merging data for all years

# setting: parameters for transforming temp to approximate normal distribution
trans_mode = 'none'  # box-cox or power-law or none
trans_exp_daily = np.nan

# infiles
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_original.npz'
file_nearstn = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_original/nearstn_catalog.npz'

# outfiles
outpath = '/home/gut428/stn_reg_original'
filereg_year = outpath + '/regression_stn_' + str(year) + '.npz'  # regression error at station points
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

    date_ymd = datatemp['date_ymd']
    yyyy = (date_ymd / 10000).astype(int)
    indy = yyyy == year
    prcp_stn_daily = datatemp['prcp_stn'][:, indy]
    tmean_stn_daily = datatemp['tmean_stn'][:, indy]
    trange_stn_daily = datatemp['trange_stn'][:, indy]

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

if not os.path.isfile(filereg_year):
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

    np.savez_compressed(filereg_year, prcp=prcp_reg_stn, tmean=tmean_reg_stn,
                        trange=trange_reg_stn, pop=pop_reg_stn)

########################################################################################################################

# check whether it is time to merge all years
flag = 1
for y in range(yearall[0], yearall[1]+1):
    filey = outpath + '/daily_regression_stn_' + str(y) + '.npz'
    if not os.path.isfile(filey):
        flag = 0
        break

if os.path.isfile(filereg_all):
    print('file for merged data exists')
else:
    print('merge data for all years')
    ndays = len(date_ymd)
    prcp = np.nan * np.zeros([nstn, ndays], dtype=np.float32)
    pop = np.nan * np.zeros([nstn, ndays], dtype=np.float32)
    tmean = np.nan * np.zeros([nstn, ndays], dtype=np.float32)
    trange = np.nan * np.zeros([nstn, ndays], dtype=np.float32)
    for y in range(yearall[0], yearall[1]+1):
        print('year',y,'/// all year', yearall)
        filey = outpath + '/daily_regression_stn_' + str(y) + '.npz'
        indy = yyyy == y
        d = np.load(filey)
        prcp[:, indy] = d['prcp']
        pop[:, indy] = d['pop']
        tmean[:, indy] = d['tmean']
        trange[:, indy] = d['trange']
        del d
    np.savez_compressed(filereg_all, prcp=prcp, tmean=tmean, trange=trange, pop=pop,
                        stninfo=stninfo, stnID=stnID, date_ymd=date_ymd)