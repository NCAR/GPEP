# interpolation from points to grids

import numpy as np
import xarray as xr
# from tqdm.contrib import itertools
import sys, os, time

from weight_calculation import distanceweight
from data_processing import data_transformation

def weighted_mean(data, weight):
    weight = weight / np.sum(weight)
    return np.dot(data, weight)


def nearby_station_statistics(stn_data, tar_nearIndex, method):
    # find every target point, find the max value of nearby estimates

    if tar_nearIndex.ndim == 2:
        # make it a 3D array to be consistent
        tar_nearIndex = tar_nearIndex[np.newaxis, :, :]

    nstn, ntime = np.shape(stn_data)
    nrow, ncol, nearmax = np.shape(tar_nearIndex)
    estimates = np.nan * np.zeros([nrow, ncol, ntime], dtype=np.float32)

    # for r, c in itertools.product(range(nrow), range(ncol)):
    for r in range(nrow):
        for c in range(ncol):
            sample_nearIndex = tar_nearIndex[r, c, :]
            index_valid = sample_nearIndex >= 0
            if np.sum(index_valid) > 0:
                sample_nearIndex = sample_nearIndex[index_valid]
                ydata_near = stn_data[sample_nearIndex, :]
                if method == 'max':
                    estimates[r, c, :] = np.nanmax(ydata_near, axis=0)
                elif method == 'min':
                    estimates[r, c, :] = np.nanmin(ydata_near, axis=0)
                else:
                    sys.exit('Unknown method!')

    return np.squeeze(estimates)


def extrapolation(datain, nearstn_loc, nearstn_DistOrWeright, weighttype, excflag=0):
    # datain: one or multiple time steps
    # excflag: excflag one value that deviates farthest from mean value. this can help reduce outliers when we have less confidence on regression
    # weighttype: 'idw', 'exp', or 'DirectWeight'. 'DirectWeight' means that nearstn_DistOrWeright is just the weight input.

    wexp = 3

    # check/change dimension
    if np.ndim(datain) == 1:  # add time axis
        datain = datain[:, np.newaxis]

    if np.ndim(nearstn_loc) == 2:
        # add a new axis to make it a 3D array
        nearstn_loc = nearstn_loc[np.newaxis, :, :]
        nearstn_DistOrWeright = nearstn_DistOrWeright[np.newaxis, :, :]

    # start interpolation
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

            if excflag == 1:
                dmean = np.tile(np.nanmean(dataini, axis=0), [np.sum(indloci), 1])
                diff = np.abs(dataini - dmean)
                for j in range(ntimes):
                    dataini[np.nanargmax(diff[:, j]), j] = np.nan

            disti = nearstn_DistOrWeright[r, c, indloci]
            if weighttype == 'exp':
                maxdist = np.max([np.max(disti) + 1, 100])
                weighti = distanceweight(disti, maxdist, wexp)
            elif weighttype == 'idw':
                disti[disti == 0] = 0.1
                weighti = 1 / (disti ** 2)
            elif weighttype == 'DirectWeight':
                weighti = disti
            else:
                sys.exit('Unknown weight type')
            weighti = weighti / np.nansum(weighti)
            weighti[np.isnan(dataini[:,0])]=np.nan
            weighti = weighti / np.nansum(weighti)
            weighti2 = np.tile(weighti, [ntimes, 1]).T
            dataout[r, c, :] = np.nansum(dataini * weighti2, axis=0)

    dataout = np.squeeze(dataout)

    return dataout


def extrapolate_auxiliary_info(config):
    t1 = time.time()

    # parse and change configurations
    case_name = config['case_name']
    outpath_parent = config['outpath_parent']
    path_regression = f'{outpath_parent}/regression/'
    os.makedirs(path_regression, exist_ok=True)

    datestamp = f"{config['date_start'].replace('-', '')}-{config['date_end'].replace('-', '')}"
    file_grid_auxiliary = f'{path_regression}/{case_name}_auxiliary_{datestamp}.nc'  # leave one out regression
    config['file_grid_auxiliary'] = file_grid_auxiliary


    # in/out information to this function
    file_allstn = config['file_allstn']
    file_cval_reg = config['file_cval_reg']
    file_stn_nearinfo = config['file_stn_nearinfo']
    file_stn_weight = config['file_stn_weight']
    file_grid_auxiliary = config['file_grid_auxiliary']

    target_vars = config['target_vars']
    if 'transform_vars' in config:
        transform_vars = config['transform_vars']
    else:
        transform_vars = [''] * len(target_vars)
        
    if 'transform' in config:
        transform_settings = config['transform']
    else:
        transform_settings = {}


    grid_lat_name = config['grid_lat_name']
    grid_lon_name = config['grid_lon_name']

    if 'target_vars_max_constrain' in config:
        target_vars_max_constrain = config['target_vars_max_constrain']
    else:
        target_vars_max_constrain = []

    if 'overwrite_stn_cv_reg' in config:
        # because error is calculated from loo
        overwrite_stn_cv_reg = config['overwrite_stn_cv_reg']
    else:
        overwrite_stn_cv_reg = False

    print('#' * 50)
    print(f'Station interpolation of CV errors')
    print('#' * 50)
    print('Input file_cval_reg:       ', file_cval_reg)
    print('Input file_stn_nearinfo:   ', file_stn_nearinfo)
    print('Input file_stn_weight:     ', file_stn_weight)
    print('Output file_grid_auxiliary:', file_grid_auxiliary)
    print('Target variables:          ', target_vars)

    if os.path.isfile(file_grid_auxiliary):
        print('Note! Output gridded error file exists')
        if overwrite_stn_cv_reg == True:
            print('overwrite_stn_cv_reg is True. Continue.')
        else:
            print('overwrite_stn_cv_reg is False. Skip regression.\n')
            return config

    ########################################################################################################################
    # initialize outputs
    with xr.open_dataset(file_cval_reg) as ds_cval:
        timeaxis = ds_cval.time.values

    with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
        xaxis = ds_nearinfo.x
        yaxis = ds_nearinfo.y

    ds_out = xr.Dataset()
    ds_out.coords['time'] = timeaxis
    ds_out.coords['x'] = xaxis
    ds_out.coords['y'] = yaxis
    ds_out[grid_lat_name] = ds_nearinfo[grid_lat_name]
    ds_out[grid_lon_name] = ds_nearinfo[grid_lon_name]

    ########################################################################################################################
    # loop variables
    for vn in range(len(target_vars)):

        var_name = target_vars[vn]
        print('interpolation of CV errors for:', var_name)

        if len(transform_vars[vn]) > 0:
            var_name_trans = var_name + '_' + transform_vars[vn]
        else:
            var_name_trans = ''

        ########################################################################################################################
        # load data

        # station data
        with xr.open_dataset(file_cval_reg) as ds_cval:
            if len(var_name_trans) > 0:
                if var_name_trans in ds_cval.data_vars:
                    loo_value = ds_cval[var_name_trans].values
                elif var_name in ds_cval.data_vars:
                    print(f'Cannot find {var_name_trans} but find {var_name} in {file_cval_reg}')
                    print(f'Apply transformation to get {var_name_trans} from {var_name}')
                    if transform_vars[vn] == 'ecdf':
                        print('ecdf is not fully supported here yet')
                    else:
                        loo_value = data_transformation(ds_cval[var_name].values, transform_vars[vn], transform_settings[transform_vars[vn]], 'transform')

                else:
                    sys.exit(f'Cannot find {var_name_trans} in {file_cval_reg}')
            else:
                loo_value = ds_cval[var_name].values

        with xr.open_dataset(file_allstn) as ds_stn:
            ds_stn = ds_stn.sel(time=ds_cval.time)
            if len(var_name_trans) > 0:
                if var_name_trans in ds_stn.data_vars:
                    stn_value = ds_stn[var_name_trans].values
                else:
                    print(f'Cannot find {var_name_trans} but find {var_name} in {file_allstn}')
                    print(f'Apply transformation to get {var_name_trans} from {var_name}')
                    if transform_vars[vn] == 'ecdf':
                        print('ecdf is not fully supported here yet')
                    else:
                        stn_value = data_transformation(ds_stn[var_name].values, transform_vars[vn], transform_settings[transform_vars[vn]], 'transform')
            else:
                stn_value = ds_stn[var_name].values

        # near information
        with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
            # nearDistance = ds_nearinfo['nearDistance_InStn_' + var_name].values
            vtmp = f'nearIndex_Grid_{var_name}'
            if vtmp in ds_nearinfo.data_vars:
                nearIndex = ds_nearinfo[vtmp].values
            else:
                sys.exit(f'Cannot find nearIndex_Grid_{var_name} in {file_stn_nearinfo}')

        # weights
        with xr.open_dataset(file_stn_weight) as ds_weight:
            vtmp = f'nearWeight_Grid_{var_name}'
            if vtmp in ds_weight.data_vars:
                nearWeight = ds_weight[vtmp].values
            else:
                sys.exit(f'Cannot find nearIndex_Grid_{var_name} in {file_stn_weight}')

        ########################################################################################################################
        # interpolation of CV errors

        # # apply a loose max/min limit to loo_value to avoid outliers in loo_value estimates
        # # this is only good when stn_value cover enough time length to be representative
        # minvalue = np.nanmin(stn_value, axis=1)[:, np.newaxis] 
        # maxvalue = np.nanmax(stn_value, axis=1)[:, np.newaxis]
        # loo_value = np.clip(loo_value, minvalue, maxvalue)  # Apply min and max constraints to loo_value to each station
        
        error = extrapolation((loo_value - stn_value) ** 2, nearIndex, nearWeight, 'DirectWeight', 0)
        error = error ** 0.5

        if len(var_name_trans) > 0:
            var_name_save = 'uncert_' + var_name_trans
        else:
            var_name_save = 'uncert_' + var_name

        ds_out[var_name_save] = xr.DataArray(error, dims=('y', 'x', 'time'))

        ########################################################################################################################
        # find the maximum observation from nearby stations
        # the maximum values will be used to calculate maximum contraints through some transformation (i.e., making it larger)
        if var_name in target_vars_max_constrain:
            print(f'Add min/max nearby values for {var_name} because it is in target_vars_max_constrain {target_vars_max_constrain}')
            estimates = nearby_station_statistics(stn_value, nearIndex, 'max')
            if len(var_name_trans) > 0:
                var_name_save = f'nearmax_{var_name_trans}'
            else:
                var_name_save = f'nearmax_{var_name}'
            if estimates.ndim == 3:
                ds_out[var_name_save] = xr.DataArray(estimates, dims=('y', 'x', 'time'))
            elif estimates.ndim == 2:
                ds_out[var_name_save] = xr.DataArray(estimates, dims=('stn', 'time'))

    # save output file
    encoding = {}
    for var in ds_out.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}

    ds_out.to_netcdf(file_grid_auxiliary, encoding=encoding)

    t2 = time.time()
    print('Time cost (seconds):', t2 - t1)
    print('Successful station interpolation of CV errors!\n\n')

    return config
