import os
import sys, time
import xarray as xr
import numpy as np
import multiprocessing

# def distance(lat1, lon1, lat2, lon2):
#     # distance from lat/lon to km
#     radius = 6371  # km
#     dlat = np.radians(lat2 - lat1)
#     dlon = np.radians(lon2 - lon1)
#     a = np.sin(dlat / 2) * np.sin(dlat / 2) + np.cos(np.radians(lat1)) \
#         * np.cos(np.radians(lat2)) * np.sin(dlon / 2) * np.sin(dlon / 2)
#     c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
#     d = radius * c
#     return d

def distance(lat1, lon1, lat2, lon2):
    # distance from lat/lon to km
    lat1r, lon1r, lat2r, lon2r = np.radians(lat1), np.radians(lon1), np.radians(lat2), np.radians(lon2)
    d = ((180 * 60) / np.pi) * (2 * np.arcsin(np.sqrt((np.sin((lat1r - lat2r) / 2)) ** 2 + np.cos(lat1r) * np.cos(lat2r) * (np.sin((lon1r - lon2r) / 2)) ** 2)))
    d = d * 1.852 # nautical mile to km
    return d

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

# parallel version
def find_nearstn_for_Grids(lat_stn, lon_stn, lat_grid, lon_grid, mask_grid, try_radius, nearstn_min, nearstn_max, initial_distance, num_processes=4):
    if lat_grid.ndim != 2:
        sys.exit('Error! Wrong dim of lat_grid!')

    # lon_stn/lat_stn can contain nan

    # simple distance threshold
    try_radius = try_radius / 100  # try within this degree (assume 1 degree ~= 100 km). if failed, expanding to all stations.

    # initialization
    nrows, ncols = np.shape(lat_grid)
    nearIndex = -99999 * np.ones([nrows, ncols, nearstn_max], dtype=int)
    nearDistance = np.nan * np.ones([nrows, ncols, nearstn_max], dtype=np.float32)

    # divide the grid into chunks for parallel processing
    chunk_size = nrows // num_processes
    chunks = [(i * chunk_size, (i + 1) * chunk_size) for i in range(num_processes)]
    chunks[-1] = (chunks[-1][0], nrows)  # last chunk may be larger if nrows is not a multiple of num_processes

    # process each chunk in parallel
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = []
        for chunk in chunks:
            result = pool.apply_async(process_chunk, (chunk, lat_stn, lon_stn, lat_grid, lon_grid, mask_grid, try_radius, nearstn_min, nearstn_max, initial_distance))
            results.append(result)

        for result, chunk in zip(results, chunks):
            chunk_ni, chunk_nd = result.get()
            nearIndex[chunk[0]:chunk[1], :, :] = chunk_ni
            nearDistance[chunk[0]:chunk[1], :, :] = chunk_nd

    return nearIndex, nearDistance

def process_chunk(chunk, lat_stn, lon_stn, lat_grid, lon_grid, mask_grid, try_radius, nearstn_min, nearstn_max, initial_distance):
    nearIndex = -99999 * np.ones([chunk[1]-chunk[0], lat_grid.shape[1], nearstn_max], dtype=int)
    nearDistance = np.nan * np.ones([chunk[1]-chunk[0], lat_grid.shape[1], nearstn_max], dtype=np.float32)

    for rr in range(chunk[0], chunk[1]):
        for cc in range(lat_grid.shape[1]):
            if mask_grid[rr, cc] == 1:
                ni, nd = find_nearstn_for_one_target(lat_grid[rr, cc], lon_grid[rr, cc], lat_stn, lon_stn, try_radius, initial_distance, nearstn_min, nearstn_max)
                nearIndex[rr - chunk[0], cc, :], nearDistance[rr - chunk[0], cc, :] = ni, nd

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


def get_near_station_info(config):

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

    num_processes = config['num_processes']

    if 'overwrite_stninfo' in config:
        overwrite_stninfo = config['overwrite_stninfo']
    else:
        overwrite_stninfo = False

    # default settings
    stn_lat_name = 'lat'
    stn_lon_name = 'lon'

    grid_lat_name = 'latitude'
    grid_lon_name = 'longitude'
    grid_mask_name = 'mask'

    print('#' * 50)
    print('Get near station information')
    print('#' * 50)
    print('input file_allstn:', file_allstn)
    print('input infile_grid_domain:', infile_grid_domain)
    print('output file_stn_nearinfo:', file_stn_nearinfo)
    print('nearstn_min:', nearstn_min)
    print('nearstn_max:', nearstn_max)
    print('try_radius:', try_radius)
    print('initial_distance:', initial_distance)
    print('Number of processes:', num_processes)

    if os.path.isfile(file_stn_nearinfo):
        print('Note! Near station info file exists')
        if overwrite_stninfo == True:
            print('overwrite_stninfo is True. Continue.')
        else:
            print('overwrite_stninfo is False. Skip finding near stations.')
            return config

    ########################################################################################################################
    # read station information

    ds_stn = xr.load_dataset(file_allstn)
    lat_stn_raw = ds_stn[stn_lat_name].values
    lon_stn_raw = ds_stn[stn_lon_name].values

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
        t11=time.time()
        nearIndex_Grid, nearDistance_Grid = find_nearstn_for_Grids(lat_stn, lon_stn, lat_grid, lon_grid, mask_grid, try_radius, nearstn_min, nearstn_max, initial_distance, num_processes)
        nearIndex_InStn, nearDistance_InStn = find_nearstn_for_InStn(lat_stn, lon_stn, try_radius, nearstn_min, nearstn_max, initial_distance)
        t22 = time.time()
        print('Time cost (seconds) of getting near station index and distance:', t22-t11)

        ########################################################################################################################
        # add near info to output file
        ds_nearinfo['nearIndex_Grid_' + vari] = xr.DataArray(nearIndex_Grid, dims=('y', 'x', 'near'))
        ds_nearinfo['nearDistance_Grid_' + vari] = xr.DataArray(nearDistance_Grid, dims=('y', 'x', 'near'))
        ds_nearinfo['nearIndex_InStn_' + vari] = xr.DataArray(nearIndex_InStn, dims=('stn', 'near'))
        ds_nearinfo['nearDistance_InStn_' + vari] = xr.DataArray(nearDistance_InStn, dims=('stn', 'near'))

    # save to output files
    encoding = {}
    for var in ds_nearinfo.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}
    ds_nearinfo.to_netcdf(file_stn_nearinfo, encoding=encoding)

    t2 = time.time()
    print('Time cost (seconds):', t2 - t1)
    print('Successful near station searching!\n\n')

    return config