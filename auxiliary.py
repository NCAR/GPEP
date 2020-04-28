import numpy as np
import sys
import netCDF4 as nc


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
    if np.ndim(destination) > 1:
        num = np.shape(destination)[0]
    else:
        num = 1
    lat1, lon1 = origin[0], origin[1]
    if num == 1:
        destination = destination[0]
        lat2, lon2 = destination[0], destination[1]
    else:
        lat2, lon2 = destination[:, 0], destination[:, 1]
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
        ind0 = (data == 0)
        datat = (data ** (1 / texp) - 1) / (1 / texp)
        datat[data == 0] = -3
    elif mode == 'power-law':
        datat = data ** (1/texp)
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
        datat = (data  / texp + 1) ** texp
    elif mode == 'power-law':
        data[data < 0] = 0
        datat = data ** texp
    elif mode == 'none':
        datat = data
    else:
        sys.exit('Unknown mode for m_transform')
    return datat


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
        near_stnWeightrc[0:nearstn_use] = (1 - (near_stnDistrc2 / max_distrc) ** 3) ** 3
        near_stnWeightrc[0:nearstn_use] = near_stnWeightrc[0:nearstn_use] / np.sum(near_stnWeightrc[0:nearstn_use])

    return near_stnLocrc, near_stnDistrc, near_stnWeightrc

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
    metout[0] = temp[0][1] # CC
    metout[1] = np.nanmean(pre - obs) # ME
    metout[2] = np.nanmean( np.abs(pre - obs)) # MAE
    metout[3] = np.sqrt(np.sum( np.square(obs - pre) ) / len(obs)) # RMSE
    return metout

