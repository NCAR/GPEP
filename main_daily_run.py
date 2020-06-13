import numpy as np
import netCDF4 as nc
import auxiliary as au
import regression as reg
import datetime as dt
from matplotlib import pyplot as plt
from scipy import io
import os
import sys

########################################################################################################################
date_cal_start = int(sys.argv[1]) # yyyymmdd
date_cal_end = int(sys.argv[2])

# 0. read/define configuration information

# setting: start and end date
# calculation start/end date:
# date_cal_start = 19790101  # yyyymmdd: start date
# date_cal_end = 19790131  # yyyymmdd: end date
# station data (in PathStn) start/end date:
date_stn_start = 19790101  # yyyymmdd: start date
date_stn_end = 20181231  # yyyymmdd: end date

# setting: paramters for lag correlation of tmean_stn_daily, and cross-correlation between prcp and trange_stn_daily
windows = 1  # parameters for auto-cc t-p-cc calculation: 1 could be better than 31
lag = 1

# setting: searching nearby stations
nearstn_min = 20  # nearby stations: minimum number
nearstn_max = 30  # nearby stations: maximum number
search_radius = 400  # km. only search stations within this radius even nearstn_max cannot be reached
max_dist = 100  # max_distance in distance-based weight calculation

# setting: parameters for transforming temp to approximate normal distribution
trans_mode = 'none'  # box-cox or power-law or none
trans_exp_daily = 4

# setting: overwrite flags. -1:don't save files; 0: don't overwrite files; 1 is to overwrite existing files;
ow_daily = 0
ow_weight = 0
ow_stn = 0

# setting: output files
datestr = str(date_cal_start) + '-' + str(date_cal_end)

## Plato settings
FileStnInfo = '/datastore/GLOBALWATER/CommonData/EMDNA/StnGridInfo/stnlist_whole.txt'  # station basic information (lists)
FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA/StnGridInfo/gridinfo_whole.nc'  # study area information
PathStn = '/home/gut428/GMET/StnInput_daily'

FileStnData = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout/stndata_' + datestr + '.npz'
FileWeight = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout/weight.npz'
FileRegError_daily = '/home/gut428/GMET/PyGMETout/error_' + datestr + '.npz'  # regression error at station points
## Plato settings

# ### Mac settings
# FileStnInfo = '/Users/localuser/GMET/pyGMET_NA/stnlist_whole.txt'  # station basic information (lists)
# FileGridInfo = '/Users/localuser/GMET/pyGMET_NA/gridinfo_whole.nc'  # study area information
# PathStn = '/Users/localuser/GMET/StnInput_daily'
#
# FileStnData = '/Users/localuser/Research/EMDNA/regression/stndata_' + datestr + '.npz'
# FileWeight = '/Users/localuser/Research/EMDNA/regression/weight_nearstn.npz'
# FileRegError_daily = '/Users/localuser/Research/EMDNA/regression/error_notrans_' + datestr + '.npz'  # regression error at station points
# FileRegError_daily_corr = '/Users/localuser/Research/EMDNA/regression/error_rescorr' + datestr + '.npz'
# FileRegression_daily = '/Users/localuser/Research/EMDNA/regression/output_' + datestr + '.npz'
# ### Mac settings

########################################################################################################################

# 1. basic information

print('Read study area basic information')
# station location and attribute information
# stninfo: [ stations, 1/lat/lon/elev/slope_ns/slope_we ]
stnID, stninfo = au.readstnlist(FileStnInfo)
nstn = len(stnID)

# time information
if date_cal_start < date_stn_start:
    sys.exit('The calculation period is earlier than the station period')
if date_cal_end > date_stn_end:
    sys.exit('The calculation period is later than the station period')

date_cal_start2 = dt.datetime.strptime(str(date_cal_start), '%Y%m%d')
date_cal_end2 = dt.datetime.strptime(str(date_cal_end), '%Y%m%d')
ntimes = (date_cal_end2 - date_cal_start2).days + 1  # time steps to be processed

date_stn_start2 = dt.datetime.strptime(str(date_stn_start), '%Y%m%d')
loc_start = (date_cal_start2 - date_stn_start2).days  # start location in the netcdf file
loc_end = loc_start + ntimes

# seconds since 1970-1-1 0:0:0
daydiff = (date_cal_start2 - dt.datetime(1970, 1, 1)).days
seconds = (np.arange(ntimes) + daydiff) * 86400

# datelist: yyyymmdd
yyyymmdd = np.zeros(ntimes, dtype=int)
for d in range(ntimes):
    dated = date_cal_start2 + dt.timedelta(days=d)
    yyyymmdd[d] = int(dated.strftime("%Y%m%d"))
yyyymm = np.floor(yyyymmdd / 100).astype(int)
mm = np.floor(np.mod(yyyymmdd, 10000) / 100).astype(int)

########################################################################################################################

# 2. read study area basic information
print('Read study area basic information')
ncfid = nc.Dataset(FileGridInfo)
gridlat = ncfid.variables['latitude'][:].data
gridlon = ncfid.variables['longitude'][:].data
gridele = ncfid.variables['elev'][:].data
gridgns = ncfid.variables['gradient_n_s'][:].data
gridgwe = ncfid.variables['gradient_w_e'][:].data
mask = ncfid.variables['mask'][:].data  # 1: grids to be considered; the other values: invalid grids
ncfid.close()

nrows, ncols = np.shape(gridlat)
gridinfo = np.zeros([nrows, ncols, 6])
gridinfo[:, :, 0] = 1
gridinfo[:, :, 1] = gridlat
gridinfo[:, :, 2] = gridlon
gridinfo[:, :, 3] = gridele
gridinfo[:, :, 4] = gridgns
gridinfo[:, :, 5] = gridgwe
del gridlat, gridlon, gridele, gridgns, gridgwe

########################################################################################################################

# 3. read data (prcp, tmin, tmax) from station files
print('Read station precipitation and temperature data')
if os.path.isfile(FileStnData) and ow_stn != 1:
    print('FileStnData exists. loading ...')
    with np.load(FileStnData) as datatemp:
        prcp_stn_daily = datatemp['prcp_stn_daily']
        tmean_stn_daily = datatemp['tmean_stn_daily']
        trange_stn_daily = datatemp['trange_stn_daily']
else:
    cai_mode = 0
    prcp_stn_daily, tmean_stn_daily, trange_stn_daily, \
    prcp_stn_climo, tmean_stn_climo, trange_stn_climo, \
    prcp_stn_anom, tmean_stn_anom, trange_stn_anom \
        = au.read_station(PathStn, stnID, loc_start, loc_end, cai_mode, yyyymm)
    np.savez_compressed(FileStnData,
                        prcp_stn_daily=prcp_stn_daily, tmean_stn_daily=tmean_stn_daily, trange_stn_daily=trange_stn_daily)
    del prcp_stn_climo, tmean_stn_climo, trange_stn_climo, prcp_stn_anom, tmean_stn_anom, trange_stn_anom

########################################################################################################################

# 5. find neighboring stations and calculate distance-based weights
if os.path.isfile(FileWeight) and ow_weight != 1:
    print('FileWeight exists. loading ...')
    with np.load(FileWeight) as datatemp:
        near_stn_prcpLoc = datatemp['near_stn_prcpLoc']
        near_stn_prcpWeight = datatemp['near_stn_prcpWeight']
        near_stn_tempLoc = datatemp['near_stn_tempLoc']
        near_stn_tempWeight = datatemp['near_stn_tempWeight']
    del datatemp

########################################################################################################################

# 6. start spatial regression

########################################################################################################################

# 6.1 estimate regression error at station points
if os.path.isfile(FileRegError_daily) and ow_daily != 1:
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

    np.savez_compressed(FileRegError_daily, pcp_reg_stn=pcp_reg_stn, tmean_reg_stn=tmean_reg_stn,
                        trange_reg_stn=trange_reg_stn, pop_reg_stn=pop_reg_stn)
