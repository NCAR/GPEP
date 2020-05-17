import numpy as np
import netCDF4 as nc
import auxiliary as au
import regression_update as reg
import datetime as dt
from matplotlib import pyplot as plt
from scipy import io
import os
import sys
import calendar


def extrapolation(datain, nearstn_loc, nearstn_weight):
    # datain: one or multiple time steps
    nrows, ncols, nearnum = np.shape(nearstn_loc)
    nstn, ntimes = np.shape(datain)
    dataout = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
    for r in range(nrows):
        for c in range(ncols):
            if not nearstn_loc[r, c, 0]>=0:
                continue
            nearloci = nearstn_loc[r, c, :]
            indloci = nearloci > -1
            dataini = datain[nearloci[indloci], :]
            weighti = nearstn_weight[r, c, indloci]
            weighti[np.isnan(dataini[:,0])]=np.nan
            weighti = weighti / np.nansum(weighti)
            weighti2 = np.tile(weighti, [ntimes, 1]).T
            dataout[r, c, :] = np.nansum(dataini * weighti2, axis=0)
    return dataout

########################################################################################################################
year = int(sys.argv[1])
for m in range(12):
    date_cal_start = year*10000 + (m + 1)*100 + 1
    date_cal_end = year*10000 + (m + 1)*100 + calendar.monthrange(year, m+1)[1]
    print('Date',date_cal_start,date_cal_end)

    # setting: output files
    datestr = str(date_cal_start) + '-' + str(date_cal_end)
    FileWeight = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout/weight.npz'
    FileRegError_daily = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout/error_' + datestr + '.npz'  # regression error at station points
    # FileRegError_daily_corr = '/home/gut428/GMET/PyGMETout/error_rescorr' + datestr + '.npz'  # regression error after residual correction
    FileRegression_daily = '/home/gut428/GMET/PyGMETout/output_realerror_' + datestr + '.npz'

    ########################################################################################################################

    # 5. find neighboring stations and calculate distance-based weights
    if os.path.isfile(FileWeight):
        print('FileWeight exists. loading ...')
        with np.load(FileWeight) as datatemp:
            near_grid_prcpLoc = datatemp['near_grid_prcpLoc']
            near_grid_prcpWeight = datatemp['near_grid_prcpWeight']
            near_grid_tempLoc = datatemp['near_grid_tempLoc']
            near_grid_tempWeight = datatemp['near_grid_tempWeight']
        del datatemp
    else:
        sys.exit('weight file does not exist')

    ########################################################################################################################

    # 6. start spatial regression

    ########################################################################################################################

    # 6.1 estimate regression error at station points
    print('FileRegError_daily exists. loading ...')
    with np.load(FileRegError_daily) as datatemp:
        pcp_err_stn_bc = datatemp['pcp_err_stn_bc']
        pcp_err_stn_raw = datatemp['pcp_err_stn_raw']
        tmean_err_stn = datatemp['tmean_err_stn']
        trange_err_stn = datatemp['trange_err_stn']
    del datatemp

    ########################################################################################################################


    # extrapolate the error to grids
    if not os.path.isfile(FileRegression_daily):
        print('estimate real errors')
        pcp_realerr_bc = extrapolation(pcp_err_stn_bc, near_grid_prcpLoc, near_grid_prcpWeight)
        pcp_realerr_raw = extrapolation(pcp_err_stn_raw, near_grid_prcpLoc, near_grid_prcpWeight)
        tmean_realerr = extrapolation(tmean_err_stn, near_grid_tempLoc, near_grid_tempWeight)
        trange_realerr = extrapolation(trange_err_stn, near_grid_tempLoc, near_grid_tempWeight)
        description = '*_realerr_* is weighted mean of estimate - observation from nearby stations. ' \
                      '*_err_* is the root mean of real errors from nearby stations'
        np.savez_compressed(FileRegression_daily, pcp_realerr_bc=pcp_realerr_bc, pcp_realerr_raw=pcp_realerr_raw,
                            trange_realerr=trange_realerr, tmean_realerr=tmean_realerr,
                            description=description)