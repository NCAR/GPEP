import os.path
import sys, time
import xarray as xr
import numpy as np

def distance(lat1, lon1, lat2, lon2):
    # distance from lat/lon to km
    radius = 6371  # km
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat / 2) * np.sin(dlat / 2) + np.cos(np.radians(lat1)) \
        * np.cos(np.radians(lat2)) * np.sin(dlon / 2) * np.sin(dlon / 2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    d = radius * c
    return d


def distanceweight(dist, maxdist = 100, exp = 3):
    if np.max(dist) > maxdist:
        maxdist = np.max(dist) + 1
    weight = (1 - (dist / maxdist) ** exp) ** exp
    return weight

def find_nearstn_for_one_target(lat_tar, lon_tar, lat_stn, lon_stn, try_radius, initial_radius, nearstn_min, nearstn_max):
    # lat_tar/lon_tar: one value
    # lat_stn/lon_stn: vector lat/lon of stations
    # try_radius: a large radius depending on the station density. this helps to reduce computation time
    # initial_radius: initial radius to find stations
    # near station number will >= nearstn_min and <= nearstn_max

    # criteria: (1) find all stations within initial_radius, (2) if the number < nearstn_min, find nearstn_min nearby
    # stations without considering initial_radius

    near_index = -99999 * np.ones(nearstn_max, dtype=int)
    near_dist = np.nan * np.ones(nearstn_max, dtype=np.float32)
    stnID = np.arange(len(lat_stn))

    # basic control to reduce the number of input stations
    try_index = (np.abs(lat_stn - lat_tar) < try_radius) & (np.abs(lon_stn - lon_tar) < try_radius)
    lat_stn_try = lat_stn[try_index]
    lon_stn_try = lon_stn[try_index]
    stnID_try = stnID[try_index]

    # calculate distance (km)
    dist_try = distance(lat_tar, lon_tar, lat_stn_try, lon_stn_try)
    index_use = (dist_try <= initial_radius)
    nstn = np.sum(index_use)
    if nstn >= nearstn_max:  # strategy-1
        dist_try = dist_try[index_use]  # delete redundant stations
        stnID_try = stnID_try[index_use]
        index_final = np.argsort(dist_try)[:nearstn_max]
        near_index[0:nearstn_max] = stnID_try[index_final]
        near_dist[0:nearstn_max] = dist_try[index_final]
    else:  # strategy-2
        dist = distance(lat_tar, lon_tar, lat_stn, lon_stn)
        index_use = dist <= initial_radius
        if np.sum(index_use) >= nearstn_min:
            stnID = stnID[index_use]
            dist = dist[index_use]
            nearstn_use = min(len(stnID), nearstn_max)
        else:
            nearstn_use = nearstn_min

        index_final = np.argsort(dist)[:nearstn_use]
        near_index[0:nearstn_use] = stnID[index_final]
        near_dist[0:nearstn_use] = dist[index_final]

    return near_index, near_dist


def find_nearstn_for_Grids(lat_stn, lon_stn, lat_grid, lon_grid, mask_grid, try_radius, nearstn_min, nearstn_max, initial_distance):
    if lat_grid.ndim != 2:
        sys.exit('Error! Wrong dim of lat_grid!')

    # lon_stn/lat_stn can contain nan

    # simple distance threshold
    try_radius = try_radius / 100  # try within this degree (assume 1 degree ~= 100 km). if failed, expanding to all stations.

    # initialization
    nrows, ncols = np.shape(lat_grid)
    nearIndex = -99999 * np.ones([nrows, ncols, nearstn_max], dtype=int)
    nearDistance = np.nan * np.ones([nrows, ncols, nearstn_max], dtype=np.float32)

    for rr in range(nrows):
        # t1 = time.time()
        # print(f'Find near station for grids. Grid row {rr} in {nrows}')
        for cc in range(ncols):
            if mask_grid[rr, cc] == 1:
                nearIndex[rr, cc, :], nearDistance[rr, cc, :] = find_nearstn_for_one_target(lat_grid[rr, cc], lon_grid[rr, cc], lat_stn, lon_stn, try_radius, initial_distance, nearstn_min, nearstn_max)
        # t2 = time.time()
        # print('Time cost (seconds):', t2-t1)

    return nearIndex, nearDistance


def find_nearstn_for_InStn(lat_stn, lon_stn, try_radius, nearstn_min, nearstn_max, initial_distance):
    # InStn: input stations themselves
    # lon_stn/lat_stn can contain nan
    # t1 = time.time()
    # print(f'Find near station for input stations')

    # simple distance threshold
    try_radius = try_radius / 100  # try within this degree (assume 1 degree ~= 100 km). if failed, expanding to all stations.

    # initialization
    nstn = len(lon_stn)
    nearIndex = -99999 * np.ones([nstn, nearstn_max], dtype=int)
    nearDistance = np.nan * np.ones([nstn, nearstn_max], dtype=np.float32)

    for i in range(nstn):
        lat_stni = lat_stn.copy()
        lon_stni = lon_stn.copy()
        lat_stni[i] = np.nan
        lon_stni[i] = np.nan
        if ~np.isnan(lat_stn[i]):
                nearIndex[i, :], nearDistance[i, :] = find_nearstn_for_one_target(lat_stn[i], lon_stn[i], lat_stni, lon_stni, try_radius, initial_distance, nearstn_min, nearstn_max)
    # t2 = time.time()
    # print('Time cost (seconds):', t2 - t1)

    return nearIndex, nearDistance


def calculate_weights_from_distance(nearDistance, initial_distance=100, exp=3):
    # calculate weights

    if nearDistance.ndim == 2:
        nstn = nearDistance.shape[0]
        nearWeight = np.nan * np.ones([nstn, nearDistance.shape[1]], dtype=np.float32)
        for i in range(nstn):
            disti = nearDistance[i, :]
            if disti[0] >= 0:
                disti = disti[disti>=0]
                max_dist = np.max([initial_distance, np.max(disti) + 1])
                nearWeight[i, 0:len(disti)] = distanceweight(disti, max_dist, exp)

    elif nearDistance.ndim == 3:
        nrow = nearDistance.shape[0]
        ncol = nearDistance.shape[1]
        nearWeight = np.nan * np.ones([nrow, ncol, nearDistance.shape[2]], dtype=np.float32)
        for i in range(nrow):
            for j in range(ncol):
                distij = nearDistance[i, j, :]
                if distij[0] >= 0:
                    distij = distij[distij>=0]
                    max_dist = np.max([initial_distance, np.max(distij) + 1])
                    nearWeight[i, j, 0:len(distij)] = distanceweight(distij, max_dist, exp)

    else:
        sys.exit('Error! nearDistance must have ndim 2 or 3.')

    return nearWeight

def get_near_station_info_and_weight(config):

    t1 = time.time()


    # parse and change configurations
    path_stn_info = config['path_stn_info']
    file_stn_nearinfo = f'{path_stn_info}/all_stn_nearinfo.nc'
    config['file_stn_nearinfo'] = file_stn_nearinfo

    # in/out information to this function
    file_allstn = config['file_allstn']
    infile_grid_domain = config['infile_grid_domain']
    file_stn_nearinfo = config['file_stn_nearinfo']

    try_radius = config['try_radius']
    initial_distance = config['initial_distance']
    nearstn_min = config['nearstn_min']
    nearstn_max = config['nearstn_max']
    # target_vars = ['prcp', 'tmean', 'trange']
    target_vars = config['target_vars']

    # default settings
    stn_lat_name = 'lat'
    stn_lon_name = 'lon'

    grid_lat_name = 'latitude'
    grid_lon_name = 'longitude'
    grid_mask_name = 'mask'


    weight_exp = 3

    print('#' * 50)
    print('Get near station information and weights')
    print('#' * 50)
    print('input file_allstn:', file_allstn)
    print('input infile_grid_domain:', infile_grid_domain)
    print('output file_stn_nearinfo:', file_stn_nearinfo)
    print('nearstn_min:', nearstn_min)
    print('nearstn_max:', nearstn_max)
    print('try_radius:', try_radius)
    print('initial_distance:', initial_distance)

    if os.path.isfile(file_stn_nearinfo):
        print('Note! Near station info file exists')
        if config['overwrite_stninfo'] == True:
            print('overwrite_stninfo is True. Continue.')
        else:
            print('overwrite_stninfo is False. Skip finding near stations.')
            return config

    ########################################################################################################################
    # read station information

    ds_stn = xr.load_dataset(file_allstn)
    lat_stn_raw = ds_stn[stn_lat_name].values
    lon_stn_raw = ds_stn[stn_lon_name].values

    # check tmean and trange. they are calculated from tmin and tmax, and thus have the same valid stations. only find near station information for tmean
    if 'tmean' in target_vars and 'trange' in target_vars:
        print('Find both tmean and trange in target_vars. Remove Trange from near station search.')
        target_vars.remove('trange')

    # for a variable, some stations do not have records
    var_mean = []
    lat_stn_valid = []
    lon_stn_valid = []
    for v in target_vars:
        vm = ds_stn[v].mean(dim='time').values
        lat_v = lat_stn_raw.copy()
        lon_v = lon_stn_raw.copy()
        lat_v[np.isnan(vm)] = np.nan
        lon_v[np.isnan(vm)] = np.nan
        var_mean.append(vm)
        lat_stn_valid.append(lat_v)
        lon_stn_valid.append(lon_v)

    ########################################################################################################################
    # read domain information
    ds_domain = xr.load_dataset(infile_grid_domain)
    # ds_domain = ds_domain.rename({'x':'lon', 'y':'lat'})
    lat_grid = ds_domain[grid_lat_name].values
    lon_grid = ds_domain[grid_lon_name].values
    mask_grid = ds_domain[grid_mask_name].values
    ds_domain.coords['y'] = lat_grid[:, 0]
    ds_domain.coords['x'] = lon_grid[0, :]

    # initialize output
    ds_nearinfo = ds_domain.copy()
    ds_nearinfo.coords['near'] = np.arange(nearstn_max)
    ds_nearinfo.coords['stn'] = ds_stn.stn.values
    for v in ds_stn.data_vars:
        if not 'time' in ds_stn[v].dims:
            ds_nearinfo['stn_'+v] = ds_stn[v]

    ########################################################################################################################
    # generate near info

    for i in range(len(target_vars)):

        lat_stn = lat_stn_valid[i]
        lon_stn = lon_stn_valid[i]
        vari = target_vars[i]

        print('Processing:', vari)

        ########################################################################################################################
        # find nearby stations for stations/grids
        nearIndex_Grid, nearDistance_Grid = find_nearstn_for_Grids(lat_stn, lon_stn, lat_grid, lon_grid, mask_grid, try_radius, nearstn_min, nearstn_max, initial_distance)
        nearIndex_InStn, nearDistance_InStn = find_nearstn_for_InStn(lat_stn, lon_stn, try_radius, nearstn_min, nearstn_max, initial_distance)

        ########################################################################################################################
        # calculate weights using distance
        nearWeight_Grid = calculate_weights_from_distance(nearDistance_Grid, initial_distance, weight_exp)
        nearWeight_InStn = calculate_weights_from_distance(nearDistance_InStn, initial_distance, weight_exp)

        ########################################################################################################################
        # add near info to output file
        ds_nearinfo['nearIndex_Grid_' + vari] = xr.DataArray(nearIndex_Grid, dims=('y', 'x', 'near'))
        ds_nearinfo['nearDistance_Grid_' + vari] = xr.DataArray(nearDistance_Grid, dims=('y', 'x', 'near'))
        ds_nearinfo['nearWeight_Grid_' + vari] = xr.DataArray(nearWeight_Grid, dims=('y', 'x', 'near'))
        ds_nearinfo['nearIndex_InStn_' + vari] = xr.DataArray(nearIndex_InStn, dims=('stn', 'near'))
        ds_nearinfo['nearDistance_InStn_' + vari] = xr.DataArray(nearDistance_InStn, dims=('stn', 'near'))
        ds_nearinfo['nearWeight_InStn_' + vari] = xr.DataArray(nearWeight_InStn, dims=('stn', 'near'))

    # save to output files
    encoding = {}
    for var in ds_nearinfo.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}
    ds_nearinfo.to_netcdf(file_stn_nearinfo, encoding=encoding)

    t2 = time.time()
    print('Time cost (seconds):', t2 - t1)
    print('Successful near station searching!')

    return config