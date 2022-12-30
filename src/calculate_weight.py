import os
import sys, time
import xarray as xr
import numpy as np

def distanceweight(dist, maxdist = 100, exp = 3):
    if np.max(dist) > maxdist:
        maxdist = np.max(dist) + 1
    weight = (1 - (dist / maxdist) ** exp) ** exp
    return weight


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

def calculate_weight_using_nearstn_info(config):

    t1 = time.time()

    # parse and change configurations
    path_stn_info = config['path_stn_info']
    file_stn_weight = f"{path_stn_info}/all_stn_weight.nc"
    config['file_stn_weight'] = file_stn_weight

    # in/out information to this function
    file_stn_nearinfo = config['file_stn_nearinfo']
    file_stn_weight = config['file_stn_weight']
    initial_distance = config['initial_distance']

    # default settings
    weight_exp = 3
    keyword = 'nearDistance' # defined in near station search
    keywords_drop = ['nearDistance', 'nearIndex']

    print('#' * 50)
    print('Calculate weights for near stations')
    print('#' * 50)
    print('input file_stn_nearinfo:', file_stn_nearinfo)
    print('output file_stn_weight:', file_stn_weight)

    if os.path.isfile(file_stn_weight):
        print('Note! Weight file exists')
        if config['overwrite_weight'] == True:
            print('overwrite_weight is True. Continue.')
        else:
            print('overwrite_weight is False. Skip weight calculation.')
            return config

    ########################################################################################################################
    # calculate weights
    ds_inout = xr.load_dataset(file_stn_nearinfo)

    for var in ds_inout.data_vars:
        if keyword in var:
            print('Processing:', var)
            nearWeight = calculate_weights_from_distance(ds_inout[var].values, initial_distance, weight_exp)
            ds_inout[var.replace(keyword, 'nearWeight')] = xr.DataArray(nearWeight, dims=ds_inout[var].dims)

    # drop some vars
    for var in ds_inout.data_vars:
        if np.any([k in var for k in keywords_drop]):
            ds_inout = ds_inout.drop_vars([var])

    # save to output files
    encoding = {}
    for var in ds_inout.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}
    ds_inout.to_netcdf(file_stn_weight, encoding=encoding)

    t2 = time.time()
    print('Time cost (seconds):', t2 - t1)
    print('Successful weight calculation!')

    return config