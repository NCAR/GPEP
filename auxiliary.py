import numpy as np
import sys
import netCDF4 as nc
import os

def distanceweight(dist,maxdist,exp):
    weight = (1 - (dist / maxdist) ** exp) ** exp
    # weight = 1 / (dist**2)
    return weight

def readstnlist(FileStnInfo):
    # read txt (station list)
    data = np.genfromtxt(FileStnInfo, delimiter=',', dtype=None, skip_header=2)
    nstn = len(data)

    stnID = [''] * nstn
    stninfo = np.zeros([nstn, 6])
    for i in range(nstn):
        stnID[i] = data[i][0].decode('utf-8')
        stninfo[i, 0] = 1  # for regression
        for j in range(5):
            stninfo[i, j + 1] = data[i][j + 1]

    return stnID, stninfo


def readstnlist_nc(FileStnInfo):
    # prototype for netcdf reading but needing modification
    ncfid = nc.Dataset(FileStnInfo)

    stnIDraw = ncfid.variables['ID'][:]
    nchar, nstn = np.shape(stnIDraw)
    # generate a list of station IDs
    stnID = [''] * nstn
    for i in range(nstn):
        temp = [''] * nchar
        for j in range(nchar):
            temp[j] = stnIDraw[j, i].decode('utf-8')
        stnID[i] = ''.join(temp)

    stninfo = np.zeros([nstn, 6])
    stninfo[:, 0] = 1  # for regression
    stninfo[:, 1] = ncfid.variables['latitude'][:].data
    stninfo[:, 2] = ncfid.variables['longitude'][:].data
    stninfo[:, 3] = ncfid.variables['elevation'][:].data
    stninfo[:, 4] = ncfid.variables['slp_ns'][:].data
    stninfo[:, 5] = ncfid.variables['slp_we'][:].data

    ncfid.close()

    return stninfo


def distance(origin, destination):
    # distance from lat/lon to km
    lat1, lon1 = origin[0], origin[1]
    if np.ndim(destination) > 1:
        lat2, lon2 = destination[:, 0], destination[:, 1]
    else:
        destination = destination[0]
        lat2, lon2 = destination[0], destination[1]

    radius = 6371  # km
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat / 2) * np.sin(dlat / 2) + np.cos(np.radians(lat1)) \
        * np.cos(np.radians(lat2)) * np.sin(dlon / 2) * np.sin(dlon / 2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    d = radius * c
    return d


def transform(data, texp, mode):
    # transform prcp to approximate normal distribution
    # mode: box-cox; power-law
    if mode == 'box-cox':
        datat = (data ** (1 / texp) - 1) / (1 / texp)
        ind0 = data == 0
        if np.sum(ind0)>0:
            datat[ind0] = -3
    elif mode == 'power-law':
        datat = data ** (1 / texp)
    elif mode == 'none':
        datat = data
    else:
        sys.exit('Unknown mode for m_transform')
    return datat


def retransform(data, texp, mode):
    # transform prcp to approximate normal distribution
    # mode: box-cox; power-law
    if mode == 'box-cox':
        data[data < -3] = -3
        datat = (data / texp + 1) ** texp
    elif mode == 'power-law':
        data[data < 0] = 0
        datat = data ** texp
    elif mode == 'none':
        datat = data
    else:
        sys.exit('Unknown mode for m_transform')
    return datat


def kge2012(obs, pre, preprocess=True):
    # Modified Kling-Gupta Efficiency (Kling et al. 2012 -
    # https://doi.org/10.1016/j.jhydrol.2012.01.011)
    # calculate error in timing and dynamics r (Pearson's correlation
    # coefficient)
    if preprocess:
        # delete the nan values
        ind_nan = np.isnan(obs) | np.isnan(pre)
        obs = obs[~ind_nan]
        pre = pre[~ind_nan]
    if len(obs) > 1:
        pre_mean = np.mean(pre, axis=0, dtype=np.float64)
        obs_mean = np.mean(obs, dtype=np.float64)
        r = np.sum((pre - pre_mean) * (obs - obs_mean), axis=0, dtype=np.float64) / \
            np.sqrt(np.sum((pre - pre_mean) ** 2, axis=0, dtype=np.float64) *
                    np.sum((obs - obs_mean) ** 2, dtype=np.float64))
        # calculate error in spread of flow gamma (avoiding cross correlation with
        # bias by dividing by the mean)
        gamma = (np.std(pre, axis=0, dtype=np.float64) / pre_mean) / \
                (np.std(obs, dtype=np.float64) / obs_mean)
        # calculate error in volume beta (bias of mean discharge)
        beta = np.mean(pre, axis=0, dtype=np.float64) / \
               np.mean(obs, axis=0, dtype=np.float64)
        # calculate the modified Kling-Gupta Efficiency KGE'
        kge = 1 - np.sqrt((r - 1) ** 2 + (gamma - 1) ** 2 + (beta - 1) ** 2)
        # kgegroup = {'KGE': kge, 'r': r, 'gamma': gamma, 'beta': beta}
        kgegroup = np.array([kge, r, gamma, beta])
    else:
        kgegroup = {'KGE': np.nan, 'r': np.nan,
                    'gamma': np.nan, 'beta': np.nan}
        kgegroup = np.array([np.nan, np.nan, np.nan, np.nan])
    return kgegroup  # or output np.vstack((KGE, r, gamma, beta))


def metric(obs, pre, preprocess=True):
    if preprocess:
        # delete the nan values
        ind_nan = np.isnan(obs) | np.isnan(pre)
        obs = obs[~ind_nan]
        pre = pre[~ind_nan]

    metout = np.zeros(4)
    temp = np.corrcoef(obs, pre)
    metout[0] = temp[0][1]  # CC
    metout[1] = np.nanmean(pre - obs)  # ME
    metout[2] = np.nanmean(np.abs(pre - obs))  # MAE
    metout[3] = np.sqrt(np.sum(np.square(obs - pre)) / len(obs))  # RMSE
    return metout


def cc_calculate(windows, lag, prcp_stn, tmean_stn, trange_stn):
    nstn, ntimes = np.shape(prcp_stn)
    if ntimes == 1:
        mean_autocorr = 0
        mean_tp_corr = 0
    else:
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
    return mean_autocorr, mean_tp_corr


def read_station(PathStn, stnID, loc_start, loc_end, cai_mode, yyyymm):
    nstn = len(stnID)
    ntimes = loc_end - loc_start

    prcp_stn_daily = np.nan * np.zeros([nstn, ntimes])
    tmin_stn_daily = np.nan * np.zeros([nstn, ntimes])
    tmax_stn_daily = np.nan * np.zeros([nstn, ntimes])
    for i in range(nstn):
        if np.mod(i,1000) == 0:
            print('Reading station data. Current station:',i,'Total stations:',nstn)
        filei = PathStn + '/' + stnID[i] + '.nc'
        try:
            ncfid = nc.Dataset(filei)
            varlist = ncfid.variables.keys()
            if 'prcp' in varlist:
                prcp_stn_daily[i, :] = ncfid.variables['prcp'][loc_start:loc_end].data
            if 'tmin' in varlist:
                tmin_stn_daily[i, :] = ncfid.variables['tmin'][loc_start:loc_end].data
            if 'tmax' in varlist:
                tmax_stn_daily[i, :] = ncfid.variables['tmax'][loc_start:loc_end].data
            ncfid.close()
        except:
            print('fail to read station:', filei)

    tmin_stn_daily[np.abs(tmin_stn_daily) > 100] = np.nan
    tmax_stn_daily[np.abs(tmax_stn_daily) > 100] = np.nan
    tmean_stn_daily = (tmin_stn_daily + tmax_stn_daily) / 2
    trange_stn_daily = np.abs(tmax_stn_daily - tmin_stn_daily)
    prcp_stn_daily[(prcp_stn_daily > 0) & (prcp_stn_daily < 0.05)] = 0  # what is the best threshold?
    del tmin_stn_daily, tmax_stn_daily

    # calculate climatological mean
    yyyymmu = np.unique(yyyymm)
    ntimes_climo = len(yyyymmu)
    prcp_stn_climo = np.zeros([nstn, ntimes_climo])
    tmean_stn_climo = np.zeros([nstn, ntimes_climo])
    trange_stn_climo = np.zeros([nstn, ntimes_climo])
    prcp_stn_anom = np.zeros([nstn, ntimes])
    tmean_stn_anom = np.zeros([nstn, ntimes])
    trange_stn_anom = np.zeros([nstn, ntimes])
    if cai_mode == 1:
        for i in range(ntimes_climo):
            indi = yyyymm == yyyymmu[i]
            prcp_stni = prcp_stn_daily[:, indi]
            tmean_stni = tmean_stn_daily[:, indi]
            trange_stni = trange_stn_daily[:, indi]

            numi = np.sum(indi).astype(int)
            # climatological mean
            if numi == 1:  # only 1 day in this month
                prcp_stn_climo[:, i] = prcp_stni
                tmean_stn_climo[:, i] = tmean_stni
                trange_stn_climo[:, i] = trange_stni

                temp = np.ones(nstn)  # equal to 1 even climo == 0
                temp[np.isnan(prcp_stni)] = np.nan
                prcp_stn_anom[:, indi] = temp

                temp = np.zeros(nstn)
                temp[np.isnan(tmean_stni)] = np.nan
                tmean_stn_anom[:, indi] = temp
                tmean_stn_anom[:, indi] = temp
            else:
                prcp_stn_climo[:, i] = np.nanmean(prcp_stni, axis=1)
                tmean_stn_climo[:, i] = np.nanmean(tmean_stni, axis=1)
                trange_stn_climo[:, i] = np.nanmean(trange_stni, axis=1)

                ind0 = prcp_stn_climo[:, i] == 0  # equal to 1 even climo == 0
                for j in range(numi):
                    temp = prcp_stni[:, j] / prcp_stn_climo[:, i]
                    temp[ind0] = 1
                    prcp_stni[:, j] = temp
                    tmean_stni[:, j] = tmean_stni[:, j] - tmean_stn_climo[:, i]
                    trange_stni[:, j] = trange_stn_daily[:, j] - trange_stn_climo[:, i]

            prcp_stn_anom[:, indi] = prcp_stni
            tmean_stn_anom[:, indi] = tmean_stni
            trange_stn_anom[:, indi] = trange_stni

    return prcp_stn_daily, tmean_stn_daily, trange_stn_daily, prcp_stn_climo, tmean_stn_climo, trange_stn_climo, \
           prcp_stn_anom, tmean_stn_anom, trange_stn_anom


def station_weight(prcp_stn_daily, tmean_stn_daily, stninfo, gridinfo, mask,
                   search_radius, nearstn_min, nearstn_max, max_dist):
    nstn = np.shape(prcp_stn_daily)[0]
    nrows, ncols, nvar = np.shape(gridinfo)
    # 5.1 for each grid cell, find its neighboring stations (near_grid_*)
    print('calculate station weights for each grid cell')

    # address stations that don't have prcp_stn_daily or tmean_stn_daily data
    latlon_prcp = np.zeros([nstn, 3])
    latlon_tmean = np.zeros([nstn, 3])
    latlon_prcp[:, 0:2] = stninfo[:, 1:3]
    latlon_tmean[:, 0:2] = stninfo[:, 1:3]
    for i in range(nstn):
        # this is for serially complete station data. if stations have missing values, this should be modified
        if np.isnan(prcp_stn_daily[i, 0]):
            latlon_prcp[i, 0:2] = np.array([np.nan, np.nan])
        if np.isnan(tmean_stn_daily[i, 0]):
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

                # prcp_stn_daily
                near_gridLocrc, near_gridDistrc, near_gridWeightrc = \
                    find_nearstn(latlon_gridrc, latlon_prcp, try_radius, search_radius, max_dist, nearstn_min,
                                 nearstn_max)
                near_grid_prcpLoc[rr, cc, :] = near_gridLocrc
                near_grid_prcpDist[rr, cc, :] = near_gridDistrc
                near_grid_prcpWeight[rr, cc, :] = near_gridWeightrc

                # tmean_stn_daily and trange_stn_daily
                near_gridLocrc, near_gridDistrc, near_gridWeightrc = \
                    find_nearstn(latlon_gridrc, latlon_tmean, try_radius, search_radius, max_dist, nearstn_min,
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
        # prcp_stn_daily
        if not np.isnan(latlon_prcp[i, 0]):
            latlon_target = latlon_prcp[i, 0:2]
            latlon_prcpi = latlon_prcp.copy()
            latlon_prcpi[i, 0:2] = np.array([np.nan, np.nan])
            near_stnLocrc, near_stnDistrc, near_stnWeightrc = \
                find_nearstn(latlon_target, latlon_prcpi, try_radius, search_radius, max_dist, nearstn_min,
                             nearstn_max)
            near_stn_prcpLoc[i, :] = near_stnLocrc
            near_stn_prcpDist[i, :] = near_stnDistrc
            near_stn_prcpWeight[i, :] = near_stnWeightrc

        # tmean_stn_daily and trange_stn_daily
        if not np.isnan(latlon_tmean[i, 0]):
            latlon_target = latlon_tmean[i, 0:2]
            latlon_tmeani = latlon_tmean.copy()
            latlon_tmeani[i, 0:2] = np.array([np.nan, np.nan])
            near_stnLocrc, near_stnDistrc, near_stnWeightrc = \
                find_nearstn(latlon_target, latlon_tmeani, try_radius, search_radius, max_dist, nearstn_min,
                             nearstn_max)
            near_stn_tempLoc[i, :] = near_stnLocrc
            near_stn_tempDist[i, :] = near_stnDistrc
            near_stn_tempWeight[i, :] = near_stnWeightrc

    return near_grid_prcpLoc, near_grid_prcpDist, near_grid_prcpWeight, \
           near_grid_tempLoc, near_grid_tempDist, near_grid_tempWeight, \
           near_stn_prcpLoc, near_stn_prcpDist, near_stn_prcpWeight, \
           near_stn_tempLoc, near_stn_tempDist, near_stn_tempWeight


def find_nearstn(latlon_gridrc, latlon_stn, try_radius, search_radius, max_dist, nearstn_min, nearstn_max):
    # criteria: (1) find all stations within search_radius, (2) if the number < nearstn_min, find nearstn_min nearby
    # stations without considering search_radius
    near_stnLocrc = -999 * np.ones(nearstn_max)
    near_stnDistrc = -999 * np.ones(nearstn_max)
    near_stnWeightrc = -999 * np.ones(nearstn_max)

    nstn = np.shape(latlon_stn)[0]
    latlon_stnrc = np.zeros([nstn, 4])
    latlon_stnrc[:, 0:3] = latlon_stn

    # basic control to reduce the number of input stations
    try_indexrc = (np.abs(latlon_stnrc[:, 0] - latlon_gridrc[0]) < try_radius) & \
                  (np.abs(latlon_stnrc[:, 1] - latlon_gridrc[1]) < try_radius)
    latlon_stnrc = latlon_stnrc[try_indexrc, :]

    # calculate distance.
    latlon_stnrc[:, 3] = distance(latlon_gridrc, latlon_stnrc[:, 0:2])

    index_use = (latlon_stnrc[:, 3] <= search_radius)
    nstnrc = np.sum(index_use)
    if nstnrc >= nearstn_max:  # strategy-1
        latlon_stnrc = latlon_stnrc[index_use, :]  # delete redundant stations
        latlon_stnrc = latlon_stnrc[np.argsort(latlon_stnrc[:, 3]), :]  # from low to high distance
        near_stnLocrc[0:nearstn_max] = latlon_stnrc[0:nearstn_max, 2]
        near_stnDistrc[0:nearstn_max] = latlon_stnrc[0:nearstn_max, 3]
    else:  # strategy-2
        latlon_stnrc = np.zeros([nstn, 4])
        latlon_stnrc[:, 0:3] = latlon_stn
        latlon_stnrc[:, 3] = distance(latlon_gridrc, latlon_stnrc[:, 0:2])
        index_use = (latlon_stnrc[:, 3] <= search_radius)
        if np.sum(index_use) >= nearstn_min:
            latlon_stnrc = latlon_stnrc[index_use, :]
            nearstn_use = min(np.shape(latlon_stnrc)[0], nearstn_max)
        else:
            nearstn_use = nearstn_min

        latlon_stnrc = latlon_stnrc[np.argsort(latlon_stnrc[:, 3]), :]  # from low to high distance
        near_stnLocrc[0:nearstn_use] = latlon_stnrc[0:nearstn_use, 2]
        near_stnDistrc[0:nearstn_use] = latlon_stnrc[0:nearstn_use, 3]

    # calculate weights
    near_stnDistrc2 = near_stnDistrc[near_stnDistrc >= 0]
    nearstn_use = np.shape(near_stnDistrc2)[0]
    if nearstn_use > 0:
        max_distrc = max(max_dist, max(near_stnDistrc2) + 1)
        temp = distanceweight(near_stnDistrc2, max_distrc, 3)
        temp[temp > max_distrc] = 0
        near_stnWeightrc[0:nearstn_use] = temp
        # near_stnWeightrc[0:nearstn_use] = near_stnWeightrc[0:nearstn_use] / np.sum(near_stnWeightrc[0:nearstn_use])

    return near_stnLocrc, near_stnDistrc, near_stnWeightrc


def save_output_nc(FileRegression, gridinfo, seconds, mean_autocorr, mean_tp_corr, pop, pcp, tmean, trange,
                   pcp_err, tmean_err, trange_err, y_max):
    # save outputs in netcdf format for probabilistic estimates
    # variables to be saved
    # pcp, pop, pcp_error, tmean, tmean_error, trange, trange_error
    # nx, ny, grdlat, grdlon, grdalt, times, mean_autocorr, mean_tp_corr, y_max

    nrows, ncols, ntimes = np.shape(pcp)

    ncfid = nc.Dataset(FileRegression, 'w', format='NETCDF4')

    ncfid.createDimension('y', nrows)
    ncfid.createDimension('x', ncols)
    ncfid.createDimension('time', ntimes)
    ncfid.createDimension('const', 1)

    varin = ncfid.createVariable('time', 'f4', ('time'), zlib=True, complevel=9)
    varin[:] = seconds
    varin.description = 'seconds since 1970-1-1 0:0:0'

    varin = ncfid.createVariable('auto_corr', 'f4', ('const'), zlib=True, complevel=9)
    varin[:] = mean_autocorr

    varin = ncfid.createVariable('tp_corr', 'f4', ('const'), zlib=True, complevel=9)
    varin[:] = mean_tp_corr

    varin = ncfid.createVariable('latitude', 'f4', ('y', 'x'), zlib=True, complevel=9)
    varin[:] = gridinfo[:, :, 1]

    varin = ncfid.createVariable('longitude', 'f4', ('y', 'x'), zlib=True, complevel=9)
    varin[:] = gridinfo[:, :, 2]

    varin = ncfid.createVariable('altitude', 'f4', ('y', 'x'), zlib=True, complevel=9)
    varin[:] = gridinfo[:, :, 3]

    varin = ncfid.createVariable('pcp', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
    dw = np.transpose(pcp, [2, 0, 1])
    varin[:] = dw

    varin = ncfid.createVariable('pop', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
    dw = np.transpose(pop, [2, 0, 1])
    varin[:] = dw

    varin = ncfid.createVariable('pcp_error', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
    dw = np.transpose(pcp_err, [2, 0, 1])
    varin[:] = dw

    varin = ncfid.createVariable('tmean', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
    dw = np.transpose(tmean, [2, 0, 1])
    varin[:] = dw

    varin = ncfid.createVariable('tmean_error', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
    dw = np.transpose(tmean_err, [2, 0, 1])
    varin[:] = dw

    varin = ncfid.createVariable('trange', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
    dw = np.transpose(trange, [2, 0, 1])
    varin[:] = dw

    varin = ncfid.createVariable('trange_error', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
    dw = np.transpose(trange_err, [2, 0, 1])
    varin[:] = dw

    varin = ncfid.createVariable('ymax', 'f4', ('time', 'y', 'x'), zlib=True, complevel=9)
    dw = np.transpose(y_max, [2, 0, 1])
    varin[:] = dw

    ncfid.close()

def climoanom_to_daily(prcp_stn_climo, tmean_stn_climo, trange_stn_climo, # station observation: climo
                       prcp_stn_anom, tmean_stn_anom, trange_stn_anom, # station observation: anomaly from the climo
                       pcp_err_stn_climo, tmean_err_stn_climo, trange_err_stn_climo, # regression error: regressed climo from station climo
                       pcp_err_stn_anom, tmean_err_stn_anom, trange_err_stn_anom, # regression error: regressed anom from station anom
                       pcp_climo, tmean_climo, trange_climo, pcp_anom, tmean_anom, trange_anom, # regression: gridded climo and anom
                       yyyymm, trans_exp_climo, trans_exp_anom, trans_mode):
    nstn, ntimes_daily = np.shape(prcp_stn_anom)
    nrows, ncols, ntimes_climo = np.shape(pcp_climo)
    yyyymmu = np.unique(yyyymm)

    prcp_stn_CAIdaily = np.zeros([nstn, ntimes_daily])
    tmean_stn_CAIdaily = np.zeros([nstn, ntimes_daily])
    trange_stn_CAIdaily = np.zeros([nstn, ntimes_daily])
    prcp_grid_CAIdaily = np.zeros([nrows, ncols, ntimes_daily])
    tmean_grid_CAIdaily = np.zeros([nrows, ncols, ntimes_daily])
    trange_grid_CAIdaily = np.zeros([nrows, ncols, ntimes_daily])

    for i in range(ntimes_climo):
        indi = yyyymm == yyyymmu[i]

        # recover at station points
        # anomaly at stations
        prcp_anomi = transform(prcp_stn_anom[:, indi], trans_exp_anom, trans_mode) + pcp_err_stn_anom[:, indi]
        prcp_anomi = retransform(prcp_anomi, trans_exp_anom, trans_mode)
        tmean_anomi = tmean_stn_anom[:, indi] + tmean_err_stn_anom[:, indi]
        trange_anomi = trange_stn_anom[:, indi] + trange_err_stn_anom[:, indi]
        # climo at stations
        prcp_climoi = transform(prcp_stn_climo[:, i],trans_exp_climo, trans_mode) + pcp_err_stn_climo[:, i]
        prcp_climoi = retransform(prcp_climoi, trans_exp_climo, trans_mode)
        tmean_climoi = tmean_stn_climo[:, i] + tmean_err_stn_climo[:, i]
        trange_climoi = trange_stn_climo[:, i] + trange_err_stn_climo[:, i]
        # combine anomaly and climo to recover daily estimates
        numi = np.sum(indi).astype(int)
        if numi == 1:
            prcp_dailyi = prcp_anomi * prcp_climoi
            tmean_dailyi = tmean_anomi + tmean_climoi
            trange_dailyi = trange_anomi + trange_climoi
        else:
            prcp_dailyi = np.zeros([nstn, numi])
            tmean_dailyi = np.zeros([nstn, numi])
            trange_dailyi = np.zeros([nstn, numi])
            for j in range(numi):
                prcp_dailyi[:, j] = prcp_anomi[:, j] * prcp_climoi
                tmean_dailyi[:, j] = tmean_anomi[:, j] + tmean_climoi
                trange_dailyi[:, j] = trange_anomi[:, j] + trange_climoi

        prcp_stn_CAIdaily[:, indi] = prcp_dailyi
        tmean_stn_CAIdaily[:, indi] = tmean_dailyi
        trange_stn_CAIdaily[:, indi] = trange_dailyi

        # recover at all grids
        prcp_anomi = pcp_anom[:, :, indi]
        prcp_anomi = retransform(prcp_anomi, trans_exp_anom, trans_mode)
        tmean_anomi = tmean_anom[:, :, indi]
        trange_anomi = trange_anom[:, :, indi]
        prcp_climoi = pcp_climo[:, :, i]
        prcp_climoi = retransform(prcp_climoi, trans_exp_climo, trans_mode)
        tmean_climoi = tmean_climo[:, :, i]
        trange_climoi = trange_climo[:, :, i]
        if numi == 1:
            prcp_dailyi = prcp_anomi * prcp_climoi
            tmean_dailyi = tmean_anomi + tmean_climoi
            trange_dailyi = trange_anomi + trange_climoi
        else:
            prcp_dailyi = np.zeros([nrows, ncols, numi])
            tmean_dailyi = np.zeros([nrows, ncols, numi])
            trange_dailyi = np.zeros([nrows, ncols, numi])
            for j in range(numi):
                prcp_dailyi[:, :, j] = prcp_anomi[:, :, j] * prcp_climoi
                tmean_dailyi[:, :, j] = tmean_anomi[:, :, j] + tmean_climoi
                trange_dailyi[:, :, j] = trange_anomi[:, :, j] + trange_climoi

        prcp_grid_CAIdaily[:, :, indi] = prcp_dailyi
        tmean_grid_CAIdaily[:, :, indi] = tmean_dailyi
        trange_grid_CAIdaily[:, :, indi] = trange_dailyi

    return prcp_stn_CAIdaily, tmean_stn_CAIdaily, trange_stn_CAIdaily, prcp_grid_CAIdaily, tmean_grid_CAIdaily, trange_grid_CAIdaily

def evaluate(pathind, prcp_stn, tmean_stn, trange_stn, # reference data (station observations)
             pcp, tmean, trange, # gridded estimates
             pcp_stnest, tmean_stnest, trange_stnest,  # estimates at station points
             stninfo, gridinfo):
    nstn, ntimes = np.shape(prcp_stn)

    # 8.1 evaluate regression for each station
    kge_stn = [0] * 3
    metric_stn = [0] * 3
    for i in range(3):
        kge_stn[i] = np.zeros([nstn, 4])
        metric_stn[i] = np.zeros([nstn, 4])

    for i in range(nstn):
        obs = prcp_stn[i, :]
        est = pcp_stnest[i, :]
        kge_stn[0][i, :] = kge2012(obs, est)
        metric_stn[0][i, :] = metric(obs, est)

        obs = tmean_stn[i, :]
        est = tmean_stnest[i, :]
        kge_stn[1][i, :] = kge2012(obs, est)
        metric_stn[1][i, :] = metric(obs, est)

        obs = trange_stn[i, :]
        est = trange_stnest[i, :]
        kge_stn[2][i, :] = kge2012(obs, est)
        metric_stn[2][i, :] = metric(obs, est)

    # 8.2 evaluate regression of gridded estimates
    kge_grid = [0] * 3
    metric_grid = [0] * 3
    for i in range(3):
        kge_grid[i] = np.zeros([nstn, 4])
        metric_grid[i] = np.zeros([nstn, 4])

    gridlat = gridinfo[:, 1, 1]
    gridlon = gridinfo[1, :, 2]
    for i in range(nstn):
        stnlat = stninfo[i, 1]
        stnlon = stninfo[i, 2]
        row = np.argmin(np.abs(stnlat - gridlat))
        col = np.argmin(np.abs(stnlon - gridlon))

        obs = prcp_stn[i, :]
        est = pcp[row, col, :]
        kge_grid[0][i, :] = kge2012(obs, est)
        metric_grid[0][i, :] = metric(obs, est)

        obs = tmean_stn[i, :]
        est = tmean[row, col, :]
        kge_grid[1][i, :] = kge2012(obs, est)
        metric_grid[1][i, :] = metric(obs, est)

        obs = trange_stn[i, :]
        est = trange[row, col, :]
        kge_grid[2][i, :] = kge2012(obs, est)
        metric_grid[2][i, :] = metric(obs, est)

    # 8.3 evaluate using independent stations
    dirtest = os.listdir(pathind)
    ntest = len(dirtest)

    kge_ind = [0] * 3
    metric_ind = [0] * 3
    for i in range(3):
        kge_ind[i] = np.zeros([ntest, 4])
        metric_ind[i] = np.zeros([ntest, 4])

    gridlat = gridinfo[:, 1, 1]
    gridlon = gridinfo[1, :, 2]
    for i in range(ntest):
        filei = pathind + '/' + dirtest[i]
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
        kge_ind[0][i, :] = kge2012(pcpobs, est)
        metric_ind[0][i, :] = metric(pcpobs, est)

        est = tmean[row, col, :]
        kge_ind[1][i, :] = kge2012(tmeanobs, est)
        metric_ind[1][i, :] = metric(tmeanobs, est)

        est = trange[row, col, :]
        kge_ind[2][i, :] = kge2012(trangeobs, est)
        metric_ind[2][i, :] = metric(trangeobs, est)

        return kge_stn, kge_grid, kge_ind, metric_stn, metric_grid, metric_ind

def stn_to_grid(data, nearloc, nearweight):
    if np.ndim(data) == 1:
        data = data[:, np.newaxis]
    nstn, ntimes = np.shape(data)
    nrows, ncols, nearnum = np.shape(nearloc)
    datagrid = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
    for r in range(nrows):
        for c in range(ncols):
            if nearloc[r, c, 0] > -1:
                locrc = nearloc[r, c, :]
                weightrc = nearweight[r, c, :]
                weightrc = weightrc[0:5]
                weightrc = weightrc / np.sum(weightrc)
                locrc = locrc[0:5]

                weightrc = np.tile(weightrc, (ntimes, 1)).T
                datarc = data[locrc, :]
                datagrid[r, c, :] = np.nansum(datarc * weightrc, axis=0)
    return datagrid