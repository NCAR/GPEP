import numpy as np
import netCDF4 as nc
import auxiliary as au
import regression as reg
import datetime as dt
from matplotlib import pyplot as plt
from scipy import io
import os
import sys
from scipy.interpolate import griddata

########################################################################################################################
# 0. read/define configuration information
# setting: file and path names of inputs
# FileStnInfo = 'stnlist_slope.w_subset.txt'  # station basic information (lists)
# FileGridInfo = 'gridinfo.0625.w_subset.nc'  # study area information
# PathStn = 'stndata'  # original station data (prcp ...)
FileStnInfo = '/Users/localuser/GMET/Example_tgq/inputs/stnlist_example.txt'  # station basic information (lists)
FileGridInfo = '/Users/localuser/GMET/Example_tgq/inputs/gridinfo_example.nc'  # study area information
# PathStn = '/Users/localuser/GMET/Example_tgq/StnDaily_train'  # original station data (prcp ...)
PathStn = '/Users/localuser/GMET/StnInput_daily'

# setting: start and end date
# calculation start/end date:
date_cal_start = 20180101  # yyyymmdd: start date
date_cal_end = 20181231  # yyyymmdd: end date
# station data (in PathStn) start/end date:
date_stn_start = 19790101  # yyyymmdd: start date
date_stn_end = 20181231  # yyyymmdd: end date

# setting: paramters for lag correlation of tmean_stn, and cross-correlation between prcp and trange_stn
windows = 31  # parameters for auto-cc t-p-cc calculation
lag = 1

# setting: searching nearby stations
nearstn_min = 10  # nearby stations: minimum number
nearstn_max = 30  # nearby stations: maximum number
search_radius = 1000  # km. only search stations within this radius even nearstn_max cannot be reached
max_dist = 100  # max_distance in distance-based weight calculation
overwrite_weight = 1  # 1: overwrite FileWeight; other values, do not overwrite FileWeight
FileWeight = '/Users/localuser/GMET/pyGMET_res/weight_nearstn.npz'

# setting: parameters for transforming temp to approximate normal distribution
trans_exp = 4
trans_mode = 'box-cox'  # box-cox or power-law or none

# setting: regression outputs
overwrite_regress = 1  # 1: overwrite  regression output files; other values, do not overwrite
FileRegError = '/Users/localuser/GMET/pyGMET_res/regress_error.npz'  # regression error at station points
FileRegression = '/Users/localuser/GMET/pyGMET_res/regress_daily.npz'

########################################################################################################################

# check file status
# this part should be activated in operational application
# if os.path.isfile(FileRegression) and overwrite_regress != 1:
#     print('Condition-1:', FileRegression, 'exists')
#     print('Condition-2: overwrite_regress != 1')
#     sys.exit('Output files have been generated. Exit the program')

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
prcp_stn = np.nan * np.zeros([nstn, ntimes])
tmin_stn = np.nan * np.zeros([nstn, ntimes])
tmax_stn = np.nan * np.zeros([nstn, ntimes])
for i in range(nstn):
    filei = PathStn + '/' + stnID[i] + '.nc'
    try:
        ncfid = nc.Dataset(filei)
        varlist = ncfid.variables.keys()
        if 'prcp' in varlist:
            prcp_stn[i, :] = ncfid.variables['prcp'][loc_start:loc_end].data
        if 'tmin' in varlist:
            tmin_stn[i, :] = ncfid.variables['tmin'][loc_start:loc_end].data
        if 'tmax' in varlist:
            tmax_stn[i, :] = ncfid.variables['tmax'][loc_start:loc_end].data
        ncfid.close()
    except:
        print('fail to read station:', filei)

tmin_stn[np.abs(tmin_stn) > 100] = np.nan
tmax_stn[np.abs(tmax_stn) > 100] = np.nan
tmean_stn = (tmin_stn + tmax_stn) / 2
trange_stn = np.abs(tmax_stn - tmin_stn)
prcp_stn[(prcp_stn > 0) & (prcp_stn < 0.05)] = 0  # what is the best threshold?
del tmin_stn, tmax_stn

########################################################################################################################

# 4. calculate auto_corr and t_p_corr
print('Calculate correlation (auto_cc and t_p_cc)')
local_window = np.min([windows, int(ntimes * 0.1)])
if lag > ntimes - 2:  # at least two values are needed to calculate cc
    print('The time lag for auto_cc calculation is too large.')
    lag = 1
    if lag > ntimes - 2:
        lag = 0

auto_corr = np.zeros(nstn)
t_p_corr = np.zeros(nstn)

for i in range(nstn):
    tmeani = tmean_stn[i, :]
    trangei = trange_stn[i, :]
    prcpi = prcp_stn[i, :]
    # smooth tmean_stn/trange_stn: moving average
    tmean_smoothi = np.convolve(tmeani, np.ones((local_window,)) / local_window, mode='same')
    trange_smoothi = np.convolve(trangei, np.ones((local_window,)) / local_window, mode='same')
    # lag auto_corr of tmean_stn
    cci = np.corrcoef(tmean_smoothi[0:-lag], tmean_smoothi[lag:])
    auto_corr[i] = cci[0, 1]
    # t_p_corr
    cci = np.corrcoef(trange_smoothi, prcpi)
    t_p_corr[i] = cci[0, 1]

auto_corr[abs(auto_corr) > 1] = np.nan
t_p_corr[abs(t_p_corr) > 1] = np.nan
mean_autocorr = np.nanmean(auto_corr)
mean_tp_corr = np.nanmean(t_p_corr)
print('Tmean lag-1 autocorrelation: ', mean_autocorr)
print('Trange-prcp correlation: ', mean_tp_corr)

########################################################################################################################

# 5. find neighboring stations and calculate distance-based weights
if os.path.isfile(FileWeight) and overwrite_weight != 1:
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
    # 5.1 for each grid cell, find its neighboring stations (near_grid_*)
    print('calculate station weights for each grid cell')
    # address stations that don't have prcp_stn or tmean_stn data
    latlon_prcp = np.zeros([nstn, 3])
    latlon_tmean = np.zeros([nstn, 3])
    latlon_prcp[:, 0:2] = stninfo[:, 1:3]
    latlon_tmean[:, 0:2] = stninfo[:, 1:3]
    for i in range(nstn):
        # this is for serially complete station data. if stations have missing values, this should be modified
        if np.isnan(prcp_stn[i, 0]):
            latlon_prcp[i, 0:2] = np.array([np.nan, np.nan])
        if np.isnan(tmean_stn[i, 0]):
            latlon_tmean[i, 0:2] = np.array([np.nan, np.nan])
    latlon_prcp[:, 2] = np.arange(nstn)  # station ID numbers/oders
    latlon_tmean[:, 2] = np.arange(nstn)  # station ID numbers/oders

    # simple distance threshold
    try_radius = search_radius / 100  # try within this degree (assume 1 degree ~= 100 km). if failed, expanding to all stations.

    # initialization
    near_grid_prcpLoc = -999 * np.ones([nrows, ncols, nearstn_max], dtype=int)
    near_grid_prcpDist = -999 * np.ones([nrows, ncols, nearstn_max], dtype=float)
    near_grid_prcpWeight = -999 * np.ones([nrows, ncols, nearstn_max], dtype=float)
    near_grid_tempLoc = -999 * np.ones([nrows, ncols, nearstn_max], dtype=int)
    near_grid_tempDist = -999 * np.ones([nrows, ncols, nearstn_max], dtype=float)
    near_grid_tempWeight = -999 * np.ones([nrows, ncols, nearstn_max], dtype=float)

    for rr in range(nrows):
        for cc in range(ncols):
            if mask[rr, cc] == 1:
                latlon_gridrc = gridinfo[rr, cc, 1:3]

                # prcp_stn
                near_gridLocrc, near_gridDistrc, near_gridWeightrc = \
                    au.find_nearstn(latlon_gridrc, latlon_prcp, try_radius, search_radius, max_dist, nearstn_min,
                                    nearstn_max)
                near_grid_prcpLoc[rr, cc, :] = near_gridLocrc
                near_grid_prcpDist[rr, cc, :] = near_gridDistrc
                near_grid_prcpWeight[rr, cc, :] = near_gridWeightrc

                # tmean_stn and trange_stn
                near_gridLocrc, near_gridDistrc, near_gridWeightrc = \
                    au.find_nearstn(latlon_gridrc, latlon_tmean, try_radius, search_radius, max_dist, nearstn_min,
                                    nearstn_max)
                near_grid_tempLoc[rr, cc, :] = near_gridLocrc
                near_grid_tempDist[rr, cc, :] = near_gridDistrc
                near_grid_tempWeight[rr, cc, :] = near_gridWeightrc

    # 5.2 for station, find its neighboring stations (near_stn_*)
    print('calculate station weights for each station')
    # initialization
    near_stn_prcpLoc = -999 * np.ones([nstn, nearstn_max], dtype=int)
    near_stn_prcpDist = -999 * np.ones([nstn, nearstn_max], dtype=float)
    near_stn_prcpWeight = -999 * np.ones([nstn, nearstn_max], dtype=float)
    near_stn_tempLoc = -999 * np.ones([nstn, nearstn_max], dtype=int)
    near_stn_tempDist = -999 * np.ones([nstn, nearstn_max], dtype=float)
    near_stn_tempWeight = -999 * np.ones([nstn, nearstn_max], dtype=float)
    for i in range(nstn):
        # prcp_stn
        if not np.isnan(latlon_prcp[i, 0]):
            latlon_target = latlon_prcp[i, 0:2]
            latlon_prcpi = latlon_prcp.copy()
            latlon_prcpi[i, 0:2] = np.array([np.nan, np.nan])
            near_stnLocrc, near_stnDistrc, near_stnWeightrc = \
                au.find_nearstn(latlon_target, latlon_prcpi, try_radius, search_radius, max_dist, nearstn_min,
                                nearstn_max)
            near_stn_prcpLoc[i, :] = near_stnLocrc
            near_stn_prcpDist[i, :] = near_stnDistrc
            near_stn_prcpWeight[i, :] = near_stnWeightrc

        # tmean_stn and trange_stn
        if not np.isnan(latlon_tmean[i, 0]):
            latlon_target = latlon_tmean[i, 0:2]
            latlon_tmeani = latlon_tmean.copy()
            latlon_tmeani[i, 0:2] = np.array([np.nan, np.nan])
            near_stnLocrc, near_stnDistrc, near_stnWeightrc = \
                au.find_nearstn(latlon_target, latlon_tmeani, try_radius, search_radius, max_dist, nearstn_min,
                                nearstn_max)
            near_stn_tempLoc[i, :] = near_stnLocrc
            near_stn_tempDist[i, :] = near_stnDistrc
            near_stn_tempWeight[i, :] = near_stnWeightrc

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
if os.path.isfile(FileRegError) and overwrite_regress != 1:
    print('FileRegError exists. loading ...')
    with np.load(FileRegError) as datatemp:
        pcp_err_stn = datatemp['pcp_err_stn']
        tmean_err_stn = datatemp['tmean_err_stn']
        trange_err_stn = datatemp['trange_err_stn']
    del datatemp
else:
    print('Estimate regression error at station points')
    pcp_err_stn = -999 * np.ones([nstn, ntimes])
    tmean_err_stn = -999 * np.ones([nstn, ntimes])
    trange_err_stn = -999 * np.ones([nstn, ntimes])

    for t in range(ntimes):
        # assign vectors of station alues for prcp_stn, temp, for current time step
        # transform prcp_stn to approximate normal distribution
        y_prcp = au.transform(prcp_stn[:, t], trans_exp, trans_mode)
        y_tmean = tmean_stn[:, t]
        y_trange = trange_stn[:, t]

        for gg in range(nstn):
            if prcp_stn[gg, t] > -1:
                # reduced matrices for precip
                nstn_prcp = int(np.sum(near_stn_prcpLoc[gg, :] > -1))
                if nstn_prcp >= nearstn_min:
                    w_pcp_red = np.zeros([nstn_prcp, nstn_prcp])
                    for i in range(nstn_prcp):
                        w_pcp_red[i, i] = near_stn_prcpWeight[gg, i]  # eye matrix: stn weight in one-one line
                    w_pcp_1d = near_stn_prcpWeight[gg, 0:nstn_prcp]  # stn weight
                    w_pcp_1d_loc = near_stn_prcpLoc[gg, 0:nstn_prcp]  # stn ID number/location
                    y_prcp_red = y_prcp[w_pcp_1d_loc]  # transformed prcp_stn
                    x_red = stninfo[w_pcp_1d_loc, :]  # station lat/lon/ele/slope_ns/slope_we

                    yp_red = np.zeros(nstn_prcp)  # pop: 0/1
                    yp_red[prcp_stn[w_pcp_1d_loc, t] > 0] = 1
                    ndata = np.sum(yp_red == 1)  # number of prcp_stn>0
                else:
                    # there are not enough nearby stations
                    x_red = 0  # not really necessary. just to stop warming from Pycharm
                    w_pcp_red = 0  # not really necessary. just to stop warming from Pycharm
                    yp_red = 0  # not really necessary. just to stop warming from Pycharm
                    y_prcp_red = 0  # not really necessary. just to stop warming from Pycharm
                    ndata = 0

                # prcp processing
                if ndata == 0:
                    # nearby stations do not have positive prcp data
                    pcp_err_stn[gg, t] = 0
                else:
                    # tmp needs to be matmul(TX, X) where TX = TWX_red and X = X_red
                    mat_test = np.matmul(np.transpose(x_red), w_pcp_red)
                    tmp = np.matmul(mat_test, x_red)
                    vv = np.max(np.abs(tmp), axis=1)

                    # decide if slope is to be used in regression
                    if np.any(vv) == 0 or (abs(stninfo[gg, 4]) < 3.6 or abs(stninfo[gg, 5]) < 3.6):
                        slope_flag_pcp = 0
                        x_red_use = x_red[:, 0:4]  # do not use slope in regression
                        stninfo_use = stninfo[gg, 0:4]
                    else:
                        slope_flag_pcp = 1
                        x_red_use = x_red
                        stninfo_use = stninfo[gg, :]

                    tx_red = np.transpose(x_red_use)
                    twx_red = np.matmul(tx_red, w_pcp_red)

                    # calculate pcp
                    b = reg.least_squares(x_red_use, y_prcp_red, twx_red)
                    pcpgg = np.dot(stninfo_use, b)
                    pcp_err_stn[gg, t] = pcpgg - y_prcp[gg]

            # tmean/trange processing
            if y_tmean[gg] > -100:
                # reduced matrices for tmean_stn/trange_stn
                nstn_temp = int(np.sum(near_stn_tempLoc[gg, :] > -1))
                if nstn_temp >= nearstn_min:
                    w_temp_red = np.zeros([nstn_temp, nstn_temp])
                    for i in range(nstn_temp):
                        w_temp_red[i, i] = near_stn_tempWeight[gg, i]  # eye matrix: stn weight in one-one lien
                    w_temp_1d = near_stn_tempWeight[gg, 0:nstn_temp]  # stn weight
                    w_temp_1d_loc = near_stn_tempLoc[gg, 0:nstn_temp]  # stn ID number/location
                    y_tmean_red = y_tmean[w_temp_1d_loc]  # transformed temp
                    y_trange_red = y_trange[w_temp_1d_loc]  # transformed temp
                    x_red_t = stninfo[w_temp_1d_loc, :]  # station lat/lon/ele/slope_ns/slope_we

                    ndata_t = np.sum(y_tmean_red > -100)
                    nodata_t = nstn_temp - ndata_t  # invalid temperature
                else:
                    x_red_t = 0  # not really necessary. just to stop warming from Pycharm
                    w_temp_red = 0  # not really necessary. just to stop warming from Pycharm
                    y_tmean_red = 0  # not really necessary. just to stop warming from Pycharm
                    y_trange_red = 0  # not really necessary. just to stop warming from Pycharm
                    ndata_t = 0
                    nodata_t = 0

                if ndata_t > 0:
                    stninfo_use = stninfo[gg, 0:4]
                    x_red_use = x_red_t[:, 0:4]  # do not use slope for temperature
                    tx_red = np.transpose(x_red_use)
                    twx_red = np.matmul(tx_red, w_temp_red)

                    b = reg.least_squares(x_red_use, y_tmean_red, tx_red)
                    tmeangg = np.dot(stninfo_use, b)
                    tmean_err_stn[gg, t] = tmeangg - y_tmean[gg]

                    b = reg.least_squares(x_red_use, y_trange_red, tx_red)
                    trangegg = np.dot(stninfo_use, b)
                    trange_err_stn[gg, t] = trangegg - y_trange[gg]

    np.savez_compressed(FileRegError, pcp_err_stn=pcp_err_stn, tmean_err_stn=tmean_err_stn,
                        trange_err_stn=trange_err_stn, stninfo=stninfo)

########################################################################################################################

# regression for each grid cell
if os.path.isfile(FileRegression) and overwrite_regress != 1:
    print('FileRegression exists. loading ...')
    with np.load(FileRegression) as datatemp:
        pop = datatemp['pop']
        pcp = datatemp['pcp']
        tmean = datatemp['tmean']
        trange = datatemp['trange']
        pcp_err = datatemp['pcp_err']
        tmean_err = datatemp['tmean_err']
        trange_err = datatemp['trange_err']
        y_max = datatemp['y_max']
    del datatemp
else:
    pass

    # 6.2 initialization
    print('Locally weighted regression of precipitation and temperature')
    tmp_weight_arr = np.eye(nearstn_max)
    y_max = -3 * np.ones([nrows, ncols, ntimes])
    pcp = -3 * np.ones([nrows, ncols, ntimes])
    pop = np.zeros([nrows, ncols, ntimes])
    pcp_err = np.zeros([nrows, ncols, ntimes])
    tmean = np.zeros([nrows, ncols, ntimes])
    trange = np.zeros([nrows, ncols, ntimes])
    tmean_err = np.zeros([nrows, ncols, ntimes])
    trange_err = np.zeros([nrows, ncols, ntimes])

    # start regression ...
    # loop through time steps
    for t in range(ntimes):
        print('Regression time step: ', t + 1, '---Total time steps: ', ntimes)
        # assign vectors of station alues for prcp_stn, temp, for current time step
        # transform prcp_stn to approximate normal distribution
        y_prcp = au.transform(prcp_stn[:, t], trans_exp, trans_mode)
        y_tmean = tmean_stn[:, t]
        y_trange = trange_stn[:, t]

        # loop through grids (row, col)
        for rr in range(nrows):
            for cc in range(ncols):
                if mask[rr, cc] != 1:
                    # the grid is outside mask extent
                    continue

                ############################################################################################################

                # 6.3 Precipitation estimation (pop and pcp)
                nstn_prcp = int(np.sum(near_grid_prcpLoc[rr, cc, :] > -1))
                if nstn_prcp < nearstn_min:
                    print('Precipitation regression: current time step, row, and col are', t, rr, cc)
                    sys.exit('Cannot find enough input stations for this grid cell')
                else:
                    # 6.3.1 reduced matrices for precipitation
                    w_pcp_red = np.zeros([nstn_prcp, nstn_prcp])
                    for i in range(nstn_prcp):
                        w_pcp_red[i, i] = near_grid_prcpWeight[rr, cc, i]  # eye matrix: stn weight in one-one lien
                    w_pcp_1d = near_grid_prcpWeight[rr, cc, 0:nstn_prcp]  # stn weight
                    w_pcp_1d_loc = near_grid_prcpLoc[rr, cc, 0:nstn_prcp]  # stn ID number/location
                    y_prcp_red = y_prcp[w_pcp_1d_loc]  # transformed prcp_stn
                    x_red = stninfo[w_pcp_1d_loc, :]  # station lat/lon/ele/slope_ns/slope_we

                    yp_red = np.zeros(nstn_prcp)  # pop: 0/1
                    yp_red[prcp_stn[w_pcp_1d_loc, t] > 0] = 1
                    ndata = np.sum(yp_red == 1)  # number of prcp_stn>0
                    nodata = np.sum(yp_red == 0)
                    y_max[rr, cc, t] = max(y_prcp_red)

                    # 6.3.2 estimate pop, pcp, pcp_err
                    # note: pcp_err is based on results from 6.1 and independent with gridded pop and pcp regression
                    if ndata == 0:
                        # nearby stations do not have positive prcp data (i.e., zero or missing)
                        pop[rr, cc, t] = 0
                        pcp[rr, cc, t] = y_prcp_red[0]  # corresponding to zero precipitation
                        pcp_err[rr, cc, t] = 0
                    else:
                        # decide if slope is to be used in regression
                        # tmp needs to be matmul(TX, X) where TX = TWX_red and X = X_red
                        mat_test = np.matmul(np.transpose(x_red), w_pcp_red)
                        tmp = np.matmul(mat_test, x_red)
                        vv = np.max(np.abs(tmp), axis=1)
                        if np.any(vv) == 0 or (abs(gridinfo[rr, cc, 4]) < 3.6 or abs(gridinfo[rr, cc, 5]) < 3.6):
                            slope_flag_pcp = 0
                            x_red_use = x_red[:, 0:4]  # do not use slope in regression
                            gridinfo_use = gridinfo[rr, cc, 0:4]
                        else:
                            slope_flag_pcp = 1
                            x_red_use = x_red
                            gridinfo_use = gridinfo[rr, cc, :]

                        tx_red = np.transpose(x_red_use)
                        twx_red = np.matmul(tx_red, w_pcp_red)

                        # calculate pop
                        if nodata == 0:
                            pop[rr, cc, t] = 1
                        else:
                            b = reg.logistic_regression(x_red_use, tx_red, yp_red)
                            zb = - np.dot(gridinfo_use, b)
                            pop[rr, cc, t] = 1 / (1 + np.exp(zb))

                        # calculate pcp
                        b = reg.least_squares(x_red_use, y_prcp_red, twx_red)
                        pcp[rr, cc, t] = np.dot(gridinfo_use, b)

                        # 6.4.3 estimate pcp error
                        err0 = pcp_err_stn[w_pcp_1d_loc, t]
                        pcp_err[rr, cc, t] = (np.sum((err0 ** 2) * w_pcp_1d) / np.sum(w_pcp_1d)) ** 0.5

                ############################################################################################################

                # 6.5 Temperature estimation (tmean and trange)
                # reduced matrices for tmean_stn/trange_stn
                nstn_temp = int(np.sum(near_grid_tempLoc[rr, cc, :] > -1))
                if nstn_temp < nearstn_min:
                    print('Temperature regression: current time step, row, and col are', t, rr, cc)
                    sys.exit('Cannot find enough input stations for this grid cell')
                else:
                    # 6.5.1 reduced matrices for precipitation
                    w_temp_red = np.zeros([nstn_temp, nstn_temp])
                    for i in range(nstn_temp):
                        w_temp_red[i, i] = near_grid_tempWeight[rr, cc, i]  # eye matrix: stn weight in one-one lien
                    w_temp_1d = near_grid_tempWeight[rr, cc, 0:nstn_temp]  # stn weight
                    w_temp_1d_loc = near_grid_tempLoc[rr, cc, 0:nstn_temp]  # stn ID number/location
                    y_tmean_red = y_tmean[w_temp_1d_loc]  # transformed temp
                    y_trange_red = y_trange[w_temp_1d_loc]  # transformed temp
                    x_red_t = stninfo[w_temp_1d_loc, :]  # station lat/lon/ele/slope_ns/slope_we

                    ndata_t = np.sum(y_tmean_red > -100)

                    if ndata_t == 0:
                        # This is not a problem for serially complete dataset (scd).
                        # But even if inputs have missing data, the way in Fortran-based GMET (simple filling) is not
                        # a good way. This problem should be solved when finding near stations.
                        print('Temperature regression: current time step, row, and col are', t, rr, cc)
                        sys.exit('Nearby stations do not have any valid temperature data')
                    else:
                        # 6.5.1 estimate tmean and its error
                        gridinfo_use = gridinfo[rr, cc, 0:4]
                        x_red_use = x_red_t[:, 0:4]  # do not use slope for temperature
                        tx_red = np.transpose(x_red_use)
                        twx_red = np.matmul(tx_red, w_temp_red)
                        b = reg.least_squares(x_red_use, y_tmean_red, tx_red)
                        tmean[rr, cc, t] = np.dot(gridinfo_use, b)

                        # error estimation
                        err0 = tmean_err_stn[w_temp_1d_loc, t]
                        tmean_err[rr, cc, t] = (np.sum((err0 ** 2) * w_pcp_1d) / np.sum(w_temp_1d)) ** 0.5

                        # 6.5.2 estimate trange and its error
                        b = reg.least_squares(x_red_use, y_trange_red, tx_red)
                        trange[rr, cc, t] = np.dot(gridinfo_use, b)

                        # error estimation
                        err0 = trange_err_stn[w_temp_1d_loc, t]
                        trange_err[rr, cc, t] = (np.sum((err0 ** 2) * w_pcp_1d) / np.sum(w_temp_1d)) ** 0.5

    np.savez_compressed(FileRegression, pop=pop, pcp=pcp, tmean=tmean, trange=trange,
                        pcp_err=pcp_err, tmean_err=tmean_err, trange_err=trange_err, y_max=y_max)
    # io.savemat('outputs.mat',{'pop': pop, 'tmean': tmean, 'trange': trange, 'pcp': pcp, 'tmean_err': tmean_err,
    #                           'trange_err': trange_err, 'pcp_err': pcp_err, 'pcp_stn': prcp_stn,
    #                           'tmean_stn': tmean_stn, 'trange_stn': trange_stn, 'stninfo': stninfo})

########################################################################################################################

# 7. save outputs
# variables to be saved
# pcp, pop, pcp_error, tmean, tmean_error, trange, trange_error
# nx, ny, grdlat, grdlon, grdalt, times, mean_autocorr, mean_tp_corr, y_max
# if (not os.path.isfile(FileRegression)) or overwrite_regress == 1:
#     ncfid = nc.Dataset(FileRegression, 'w', format='NETCDF4')
#
#     ncfid.createDimension('y', nrows)
#     ncfid.createDimension('x', ncols)
#     ncfid.createDimension('time', ntimes)
#     ncfid.createDimension('const', 1)
#
#     varin = ncfid.createVariable('time', 'f4', ('time'), zlib=True, complevel=9)
#     varin[:] = seconds
#     varin.description = 'seconds since 1970-1-1 0:0:0'
#
#     varin = ncfid.createVariable('auto_corr', 'f4', ('const'), zlib=True, complevel=9)
#     varin[:] = mean_autocorr
#
#     varin = ncfid.createVariable('tp_corr', 'f4', ('const'), zlib=True, complevel=9)
#     varin[:] = mean_tp_corr
#
#     varin = ncfid.createVariable('latitude', 'f4', ('y', 'x'), zlib=True, complevel=9)
#     varin[:] = gridinfo[:, :, 1]
#
#     varin = ncfid.createVariable('longitude', 'f4', ('y', 'x'), zlib=True, complevel=9)
#     varin[:] = gridinfo[:, :, 2]
#
#     varin = ncfid.createVariable('altitude', 'f4', ('y', 'x'), zlib=True, complevel=9)
#     varin[:] = gridinfo[:, :, 3]
#
#     varin = ncfid.createVariable('pcp', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
#     dw = np.transpose(pcp, [2, 0, 1])
#     varin[:] = dw
#
#     varin = ncfid.createVariable('pop', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
#     dw = np.transpose(pop, [2, 0, 1])
#     varin[:] = dw
#
#     varin = ncfid.createVariable('pcp_error', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
#     dw = np.transpose(pcp_err, [2, 0, 1])
#     varin[:] = dw
#
#     varin = ncfid.createVariable('tmean', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
#     dw = np.transpose(tmean, [2, 0, 1])
#     varin[:] = dw
#
#     varin = ncfid.createVariable('tmean_error', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
#     dw = np.transpose(tmean_err, [2, 0, 1])
#     varin[:] = dw
#
#     varin = ncfid.createVariable('trange', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
#     dw = np.transpose(trange, [2, 0, 1])
#     varin[:] = dw
#
#     varin = ncfid.createVariable('trange_error', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
#     dw = np.transpose(trange_err, [2, 0, 1])
#     varin[:] = dw
#
#     varin = ncfid.createVariable('ymax', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
#     dw = np.transpose(y_max, [2, 0, 1])
#     varin[:] = dw
#
#     ncfid.close()

########################################################################################################################

# 8. Evaluate regression results
# this part should be moved for the operational version
# 8.1 evaluate regression for each station
kge_stn = [0] * 3
metric_stn = [0]*3
kge_pcp_t = [0] * 3 # transformed pcp
kge_pcp_t[0] = np.zeros([nstn, 4])
for i in range(3):
    kge_stn[i] = np.zeros([nstn, 4])
    metric_stn[i] = np.zeros([nstn, 4])

for i in range(nstn):
    obs = prcp_stn[i, :].copy()
    obst = au.transform(obs, trans_exp, trans_mode)
    est = au.retransform(obst + pcp_err_stn[i, :], trans_exp, trans_mode)
    kge_stn[0][i, :] = au.kge2012(obs, est)
    metric_stn[0][i, :] = au.metric(obs, est)

    obs = prcp_stn[i, :].copy()
    obst = au.transform(obs, trans_exp, trans_mode)
    est = obst + pcp_err_stn[i, :]
    kge_pcp_t[0][i, :] = au.kge2012(obst, est)


    obs = tmean_stn[i, :]
    est = obs + tmean_err_stn[i, :]
    kge_stn[1][i, :] = au.kge2012(obs, est)
    metric_stn[1][i, :] = au.metric(obs, est)

    obs = trange_stn[i, :]
    est = obs + trange_err_stn[i, :]
    kge_stn[2][i, :] = au.kge2012(obs, est)
    metric_stn[2][i, :] = au.metric(obs, est)

# 8.2 evaluate regression of gridded estimates
kge_grid = [0] * 3
for i in range(3):
    kge_grid[i] = np.zeros([nstn, 4])
kge_pcp_t[1] = np.zeros([nstn, 4])

gridlat = gridinfo[:, 1, 1]
gridlon = gridinfo[1, :, 2]
for i in range(nstn):
    stnlat = stninfo[i, 1]
    stnlon = stninfo[i, 2]
    row = np.argmin(np.abs(stnlat - gridlat))
    col = np.argmin(np.abs(stnlon - gridlon))

    obs = prcp_stn[i, :].copy()
    est = pcp[row, col, :].copy()
    est = au.retransform(est, trans_exp, trans_mode)
    kge_grid[0][i, :] = au.kge2012(obs, est)

    obs = prcp_stn[i, :].copy()
    obs = au.transform(obs, trans_exp, trans_mode)
    est = pcp[row, col, :].copy()
    kge_pcp_t[1][i, :] = au.kge2012(obs, est)

    obs = tmean_stn[i, :]
    est = tmean[row, col, :]
    kge_grid[1][i, :] = au.kge2012(obs, est)

    obs = trange_stn[i, :]
    est = trange[row, col, :]
    kge_grid[2][i, :] = au.kge2012(obs, est)

# 8.3 evaluate using independent stations
pathtest = '/Users/localuser/GMET/Example_tgq/StnDaily_test'
dirtest = os.listdir(pathtest)
ntest = len(dirtest)

kge_ind = [0] * 3
for i in range(3):
    kge_ind[i] = np.zeros([ntest, 4])
kge_pcp_t[2] = np.zeros([ntest, 4])

gridlat = gridinfo[:, 1, 1]
gridlon = gridinfo[1, :, 2]
for i in range(ntest):
    filei = pathtest + '/' + dirtest[i]
    ncfid = nc.Dataset(filei)
    stnlat = ncfid.variables['latitude'][:]
    stnlon = ncfid.variables['longitude'][:]

    varlist = ncfid.variables.keys()
    if 'prcp' in varlist:
        pcpobs = ncfid.variables['prcp'][:].data
    else:
        pcpobs = np.nan * np.ones(ntimes)
    if 'tmin' in varlist:
        tminobs = ncfid.variables['tmin'][:].data
    else:
        tminobs = np.nan * np.ones(ntimes)
    if 'tmax' in varlist:
        tmaxobs = ncfid.variables['tmax'][:].data
    else:
        tmaxobs = np.nan * np.ones(ntimes)
    ncfid.close()

    row = np.argmin(np.abs(stnlat - gridlat))
    col = np.argmin(np.abs(stnlon - gridlon))
    tmeanobs = (tminobs + tmaxobs) / 2
    trangeobs = (tmaxobs - tminobs)

    est = pcp[row, col, :]
    est = au.retransform(est, trans_exp, trans_mode)
    kge_ind[0][i, :] = au.kge2012(pcpobs, est)

    pcpobs = au.transform(pcpobs, trans_exp, trans_mode)
    est = pcp[row, col, :]
    kge_pcp_t[2][i, :] = au.kge2012(pcpobs, est)

    est = tmean[row, col, :]
    kge_ind[1][i, :] = au.kge2012(tmeanobs, est)

    est = trange[row, col, :]
    kge_ind[2][i, :] = au.kge2012(trangeobs, est)


print('median kge')
print('-------------------------------------')
print('pcp after transformation')
print('station, grid_train, grid_test' )
print('%.3f  %.3f   %.3f' % (np.nanmedian(kge_pcp_t[0][:,0]), np.nanmedian(kge_pcp_t[1][:,0]), np.nanmedian(kge_pcp_t[2][:,0])))
print('-------------------------------------')
print('pcp before transformation')
print('station, grid_train, grid_test' )
print('%.3f  %.3f   %.3f' % (np.nanmedian(kge_stn[0][:,0]), np.nanmedian(kge_grid[0][:,0]), np.nanmedian(kge_ind[0][:,0])))
print('-------------------------------------')
print('tmean')
print('station, gridtrain, gridtest' )
print('%.3f  %.3f   %.3f' % (np.nanmedian(kge_stn[1][:,0]), np.nanmedian(kge_grid[1][:,0]), np.nanmedian(kge_ind[1][:,0])))
print('-------------------------------------')
print('trange')
print('tation, grid_train, grid_test' )
print('%.3f  %.3f   %.3f' % (np.nanmedian(kge_stn[2][:,0]), np.nanmedian(kge_grid[2][:,0]), np.nanmedian(kge_ind[2][:,0])))