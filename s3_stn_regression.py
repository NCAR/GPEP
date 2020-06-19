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
date_cal_start = year*10000+100+1
date_cal_end = year*10000+1200+31

# 0. read/define configuration information

# setting: start and end date
# calculation start/end date:
# date_cal_start = 19790101  # yyyymmdd: start date
# date_cal_end = 19790131  # yyyymmdd: end date
# station data (in PathStn) start/end date:
date_stn_start = 19790101  # yyyymmdd: start date
date_stn_end = 20181231  # yyyymmdd: end date

# setting: searching nearby stations
nearstn_min = 20  # nearby stations: minimum number
nearstn_max = 30  # nearby stations: maximum number
search_radius = 400  # km. only search stations within this radius even nearstn_max cannot be reached
max_dist = 100  # max_distance in distance-based weight calculation

# setting: parameters for transforming temp to approximate normal distribution
trans_mode = 'none'  # box-cox or power-law or none
trans_exp_daily = 4

# setting: output files
datestr = str(date_cal_start) + '-' + str(date_cal_end)

## Plato settings
FileStnInfo = '/datastore/GLOBALWATER/CommonData/EMDNA/StnGridInfo/stnlist_whole.txt'  # station basic information (lists)
FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA/StnGridInfo/gridinfo_whole.nc'  # study area information
PathStn = '/home/gut428/GMET/StnInput_daily'

gmet_stndatafile = '/home/gut428/stndata_whole.npz'
FileWeight = '/home/gut428/GMET/PyGMETout/weight.npz'
FileRegError_daily = '/home/gut428/GMET/PyGMETout/error_' + datestr + '.npz'  # regression error at station points
## Plato settings

########################################################################################################################

# load station data
if not os.path.isfile(gmet_stndatafile):
    print(gmet_stndatafile,'does not exist')
    sys.exit()
else:
    print('load station data')
    datatemp=np.load(gmet_stndatafile)
    date_number = datatemp['date_number']
    indy = date_number['yyyy'] == year
    prcp_stn_daily = datatemp['prcp_stn'][:, indy]
    tmean_stn_daily = datatemp['tmean_stn'][:, indy]
    trange_stn_daily = datatemp['trange_stn'][:, indy]

########################################################################################################################

# find neighboring stations and calculate distance-based weights
if os.path.isfile(FileWeight):
    print('FileWeight exists. loading ...')
    with np.load(FileWeight) as datatemp:
        near_grid_prcpLoc = datatemp['near_grid_prcpLoc']
        near_grid_prcpWeight = datatemp['near_grid_prcpWeight']
        near_grid_tempLoc = datatemp['near_grid_tempLoc']
        near_grid_tempWeight = datatemp['near_grid_tempWeight']
        near_stn_prcpLoc = datatemp['near_stn_prcpLoc']
        near_stn_prcpWeight = datatemp['near_stn_prcpWeight']
        near_stn_tempLoc = datatemp['near_stn_tempLoc']
        near_stn_tempWeight = datatemp['near_stn_tempWeight']
    del datatemp
else:
    near_grid_prcpLoc, near_grid_prcpDist, near_grid_prcpWeight, \
    near_grid_tempLoc, near_grid_tempDist, near_grid_tempWeight, \
    near_stn_prcpLoc, near_stn_prcpDist, near_stn_prcpWeight, \
    near_stn_tempLoc, near_stn_tempDist, near_stn_tempWeight \
        = au.station_weight(prcp_stn_daily, tmean_stn_daily, stninfo, gridinfo, mask,
                            search_radius, nearstn_min, nearstn_max, max_dist)

    # save data
    np.savez_compressed(FileWeight, near_grid_prcpLoc=near_grid_prcpLoc, near_grid_prcpDist=near_grid_prcpDist,
                        near_grid_prcpWeight=near_grid_prcpWeight, near_grid_tempLoc=near_grid_tempLoc,
                        near_grid_tempDist=near_grid_tempDist, near_grid_tempWeight=near_grid_tempWeight,
                        near_stn_prcpLoc=near_stn_prcpLoc, near_stn_prcpDist=near_stn_prcpDist,
                        near_stn_prcpWeight=near_stn_prcpWeight, near_stn_tempLoc=near_stn_tempLoc,
                        near_stn_tempDist=near_stn_tempDist, near_stn_tempWeight=near_stn_tempWeight)

########################################################################################################################

# 6. start spatial regression

########################################################################################################################

# 6.1 estimate regression error at station points
if os.path.isfile(FileRegError_daily):
    print('FileRegError_daily exists. loading ...')
    with np.load(FileRegError_daily) as datatemp:
        pcp_err_stn_daily = datatemp['pcp_err_stn']
        tmean_err_stn_daily = datatemp['tmean_err_stn']
        trange_err_stn_daily = datatemp['trange_err_stn']
    del datatemp
else:
    print('Estimate daily regression error at station points')
    pcp_err_stn_daily, tmean_err_stn_daily, trange_err_stn_daily, pop_err_stn_daily = \
        reg.station_error(prcp_stn_daily, tmean_stn_daily, trange_stn_daily, stninfo, near_stn_prcpLoc,
                          near_stn_prcpWeight, near_stn_tempLoc, near_stn_tempWeight, trans_exp_daily,
                          trans_mode, nearstn_min)

    pcp_reg_stn = pcp_err_stn_daily + prcp_stn_daily
    tmean_reg_stn = tmean_err_stn_daily + tmean_stn_daily
    trange_reg_stn = trange_err_stn_daily + trange_stn_daily

    prcp_stn_daily[prcp_stn_daily>0] = 1
    pop_reg_stn = pop_err_stn_daily + prcp_stn_daily

    np.savez_compressed(FileRegError_daily, pcp=pcp_reg_stn, tmean=tmean_reg_stn,
                        trange=trange_reg_stn, pop=pop_reg_stn)