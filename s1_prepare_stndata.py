# prepare station data for later use and perform preliminary quality control
import numpy as np
import os
import netCDF4 as nc
import datetime
import pandas as pd
import auxiliary as au

########################################################################################################################
# basic settings

# input files that meet GMET input formats
gmet_stnpath = '/Users/guoqiang/Github/GMET/test_cases/cali2017/stndata'  # station files for gmet_stnfile
gmet_stnfile = '/Users/guoqiang/Github/GMET/test_cases/cali2017/inputs/CALI.screened_stn_list_slope.v3.csv'  # station lists

date_start = datetime.date(2017, 2, 1)
date_end = datetime.date(2017, 2, 28)

minTrange = 1 # minimum temperature range allowed

# output file
outpath = './Andy_test_case'
os.makedirs(outpath, exist_ok=True)
gmet_stndatafile = f'{outpath}/stndata_whole.npz'  # to be saved. only process when absent

########################################################################################################################
# start processing

if os.path.isfile(gmet_stndatafile):
    print('File already exists')
else:
    # generate date list
    daynum = (date_end - date_start).days + 1
    date_ymd = np.zeros(daynum, dtype=int)
    dated = date_start
    for d in range(daynum):
        if d > 0:
            dated = dated + datetime.timedelta(days=1)
        date_ymd[d] = int(dated.strftime("%Y%m%d"))

    # read station information
    stnID, stninfo = au.readstnlist(gmet_stnfile)
    nstn = len(stnID)

    # read data from input path according to station list
    prcp_stn = np.nan * np.zeros([nstn, daynum], dtype=np.float32)
    tmin_stn = np.nan * np.zeros([nstn, daynum], dtype=np.float32)
    tmax_stn = np.nan * np.zeros([nstn, daynum], dtype=np.float32)
    for i in range(nstn):
        if np.mod(i, 1000) == 0:
            print('reading station', i, nstn)
        file = gmet_stnpath + '/' + stnID[i] + '.nc'
        fid = nc.Dataset(file)
        varlist = fid.variables.keys()
        if 'prcp' in varlist:
            prcp_stn[i, :] = fid['prcp'][:].data
        if 'tmin' in varlist:
            tmin_stn[i, :] = fid['tmin'][:].data
        if 'tmax' in varlist:
            tmax_stn[i, :] = fid['tmax'][:].data
        fid.close()
    tmean_stn = (tmin_stn + tmax_stn) / 2
    trange_stn = np.abs(tmax_stn - tmin_stn)
    trange_stn[trange_stn < minTrange] = minTrange

    # basic control of prcp_stn (cannot be too dry or just has very few unique values)
    for i in range(nstn):
        if np.isnan(prcp_stn[i, 0]):
            continue
        if (np.sum(prcp_stn[i, :]>0) / np.shape(prcp_stn)[1]) < 0.005:
            prcp_stn[i, :] = np.nan
        if len(np.unique(prcp_stn[i, :])) < 15:
            prcp_stn[i, :] = np.nan

    # save output data
    np.savez_compressed(gmet_stndatafile, prcp_stn=prcp_stn, tmean_stn=tmean_stn, trange_stn=trange_stn,
                        date_ymd=date_ymd, stnID=stnID, stninfo=stninfo)