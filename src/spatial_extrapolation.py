# interpolation from points to grids

import numpy as np
import xarray as xr
import sys, os

from weight_calculation import distanceweight

def weighted_mean(data, weight):
    weight = weight / np.sum(weight)
    return np.dot(data, weight)



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


def station_error_extrapolation(config):

    # parse and change configurations
    outpath_parent = config['outpath_parent']
    path_regression = f'{outpath_parent}/regression_outputs'
    os.makedirs(path_regression, exist_ok=True)
    file_grid_error = f'{path_regression}/Grid_Error.nc'  # leave one out regression
    config['file_grid_error'] = file_grid_error


    # in/out information to this function
    file_allstn = config['file_allstn']
    file_loo_reg = config['file_loo_reg']
    file_stn_nearinfo = config['file_stn_nearinfo']
    file_stn_weight = config['file_stn_weight']
    file_grid_error = config['file_grid_error']

    target_vars = config['target_vars']
    transform_vars = config['transform_vars']
    overwrite_loo_reg = config['overwrite_loo_reg']

    print('#' * 50)
    print(f'Station error interpolation')
    print('#' * 50)
    print('Input file_loo_reg:', file_loo_reg)
    print('Input file_stn_nearinfo:', file_stn_nearinfo)
    print('Input file_stn_weight:', file_stn_weight)
    print('Output file_grid_error:', file_grid_error)
    print('Target variables:', target_vars)

    if os.path.isfile(file_grid_error):
        print('Note! Output gridded error file exists')
        if overwrite_loo_reg == True:
            print('overwrite_loo_reg is True. Continue.')
        else:
            print('overwrite_loo_reg is False. Skip regression.')
            return config

    ########################################################################################################################
    # initialize outputs
    with xr.open_dataset(file_loo_reg) as ds_loo:
        timeaxis = ds_loo.time.values

    with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
        xaxis = ds_nearinfo.x
        yaxis = ds_nearinfo.y

    ds_out = xr.Dataset()
    ds_out.coords['time'] = timeaxis
    ds_out.coords['x'] = xaxis
    ds_out.coords['y'] = yaxis

    ########################################################################################################################
    # loop variables
    for vn in range(len(target_vars)):

        var_name = target_vars[vn]
        print('Error interpolation for:', var_name)

        if len(transform_vars[vn]) > 0:
            var_name_trans = var_name + '_' + transform_vars[vn]
        else:
            var_name_trans = ''

        ########################################################################################################################
        # load data

        # station data
        with xr.open_dataset(file_loo_reg) as ds_loo:
            if len(var_name_trans) > 0:
                loo_value = ds_loo[var_name_trans].values
            else:
                loo_value = ds_loo[var_name].values

        with xr.open_dataset(file_allstn) as ds_stn:
            ds_stn = ds_stn.sel(time=ds_loo.time)
            if len(var_name_trans) > 0:
                stn_value = ds_stn[var_name_trans].values
            else:
                stn_value = ds_stn[var_name].values

        # near information
        with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
            # nearDistance = ds_nearinfo['nearDistance_InStn_' + var_name].values
            vtmp = f'nearIndex_Grid_{var_name}'
            if vtmp in ds_nearinfo.data_vars:
                nearIndex = ds_nearinfo[vtmp].values
            else:
                if var_name == 'trange':
                    print(f'Use nearIndex_Grid_tmean for trange')
                    nearIndex = ds_nearinfo[f'nearIndex_Grid_tmean'].values
                else:
                    sys.exit(f'Cannot find nearIndex_Grid_{var_name} in {file_stn_nearinfo}')

        # weights
        with xr.open_dataset(file_stn_weight) as ds_weight:
            vtmp = f'nearWeight_Grid_{var_name}'
            if vtmp in ds_weight.data_vars:
                nearWeight = ds_weight[vtmp].values
            else:
                if var_name == 'trange':
                    print(f'Use nearWeight_Grid_tmean for trange')
                    nearWeight = ds_weight[f'nearWeight_Grid_tmean'].values
                else:
                    sys.exit(f'Cannot find nearIndex_Grid_{var_name} in {file_stn_weight}')

        ########################################################################################################################
        # error interpolation
        error = extrapolation((loo_value - stn_value) ** 2, nearIndex, nearWeight, 'DirectWeight', 0)
        error = error ** 0.5

        if len(var_name_trans) > 0:
            var_name_save = var_name_trans
        else:
            var_name_save = var_name

        ds_out[var_name_save] = xr.DataArray(error, dims=('y', 'x', 'time'))

    # save output file
    encoding = {}
    for var in ds_out.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}

    ds_out.to_netcdf(file_grid_error, encoding=encoding)

    return config
