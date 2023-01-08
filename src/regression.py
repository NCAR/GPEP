import time

import numpy as np
import xarray as xr
from tqdm.contrib import itertools
from data_processing import data_transformation
import spatial_extrapolation
import sys, os

########################################################################################################################
# basic regression utility functions

def least_squares(x, y, tx):
    # In fortran version, ludcmp and lubksb are used to calcualte matrix inversion
    # call ludcmp(a, indx, d)
    # call lubksb(a, indx, b)

    # In Python version, numpy is used to calculate matrix inversion
    c = np.matmul(tx, y)
    a = np.matmul(tx, x)

    n = np.shape(a)[0]
    b = np.zeros(n)

    deta = np.linalg.det(a)  # Compute the determinant of an array
    if deta == 0:
        # print('Singular matrix')
        b[:] = 0
    else:
        ainv = np.linalg.inv(a)
        b = np.matmul(ainv, c)

    return b

def logistic_regression(x, tx, yp):
    # nstn: station number
    # nvars: station attributes (1, lat/lon/...), 1 is for regression
    nstn, nvars = np.shape(x)

    b = np.zeros(nvars)  # regression coefficients (beta)
    p = np.zeros(nstn)  # estimated pop (probability of precipitation)
    f = 0  # flag: continue or stop loops
    it = 0  # iteration times

    while f != 1:
        # check for divergence
        xb = -np.matmul(x, b)
        if np.any(xb > 50):
            f = 1
        else:
            p = 1 / (1 + np.exp(xb))

        # check for divergence
        if np.any(p > 0.9999):
            # logistic regression diverging
            f = 1
        else:
            v = np.zeros([nstn, nstn])  # diagonal variance matrix
            for i in range(nstn):
                v[i, i] = p[i] * (1 - p[i])
            xv = np.matmul(v, x)
            yn = yp - p  # difference between station occurrence (0/1) and estimated pop: Bnew -Bold in Martyn and Slater 2006
            bn = least_squares(xv, yn, tx)

            # check: converging
            if np.any(np.abs(bn) > 1e-4):
                f = 0
            else:
                f = 1

            # check: iteration times
            if it > 20:
                f = 1
            else:
                f = 0

            b = b + bn  # update coefficients

        it = it + 1

    return b


def weight_linear_regression(nearinfo, weightnear, datanear, tarinfo):
    # # nearinfo: predictors from neighboring stations
    # # [station number, predictor number + 1] array with the first column being ones
    # nearinfo = np.zeros([nnum, npred+1])
    #
    # # weightnear: weight of neighboring stations
    # # [station number, station number] array with weights located in the diagonal
    # weightnear = np.zeros([nnum, nnum])
    # for i in range(nnum):
    #     weightnear[i, i] = 123
    #
    # # tarinfo:  predictors from target stations
    # # [predictor number + 1] vector with the first value being one
    # tarinfo = np.zeros(npred+1)
    #
    # # datanear: data from neighboring stations. [station number] vector
    # datanear = np.zeros(nnum)

    # start regression
    w_pcp_red = np.diag(np.squeeze(weightnear))
    tx_red = np.transpose(nearinfo)
    twx_red = np.matmul(tx_red, w_pcp_red)
    b = least_squares(nearinfo, datanear, twx_red)
    datatar = np.dot(tarinfo, b)

    return datatar

def weight_logistic_regression(nearinfo, weightnear, datanear, tarinfo):

    w_pcp_red = np.diag(np.squeeze(weightnear))

    # if len(np.unique(datanear)) == 1:
    #     pop = datanear[0] # all 0 (no rain) or all 1 (rain everywhere)
    # else:
    tx_red = np.transpose(nearinfo)
    twx_red = np.matmul(tx_red, w_pcp_red)
    b = logistic_regression(nearinfo, twx_red, datanear)
    if np.all(b == 0) or np.any(np.isnan(b)):
        pop = np.dot(weightnear, datanear)
    else:
        zb = - np.dot(tarinfo, b)
        pop = 1 / (1 + np.exp(zb))

    return pop

###########################
# regression using Python sklearn package: easier to use, but slower and a bit different from the above functions
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
def sklearn_weight_linear_regression(nearinfo, weightnear, datanear, tarinfo):
    model = LinearRegression()
    model = model.fit(nearinfo, datanear, sample_weight=weightnear)
    datatar = model.predict(tarinfo[np.newaxis, :])
    return datatar
def sklearn_weight_logistic_regression(nearinfo, weightnear, datanear, tarinfo):

    model = LogisticRegression(solver='liblinear')
    model = model.fit(nearinfo, datanear, sample_weight=weightnear)
    pop = model.predict_proba(tarinfo[np.newaxis,:])[0, 1]
    return pop


########################################################################################################################
# wrapped up regression functions

def loop_regression_2Dor3D(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight, tar_predictor, method):
    # regression for 2D (vector) or 3D (array) inputs
    # 2D:
    # stn_data: [input station, time steps]
    # stn_predictor: [input station, number of predictors]
    # nearIndex/nearWeight: [target point, number of nearby stations]
    # tar_predictor: [target point, number of predictors]
    # 3D:
    # stn_data: [input station, time steps]
    # stn_predictor: [input station, number of predictors]
    # nearIndex/nearWeight: [row, col, number of nearby stations]
    # tar_predictor: [row, col, number of predictors]

    if tar_nearIndex.ndim == 2:
        # make it a 3D array to be consistent
        tar_nearIndex = tar_nearIndex[np.newaxis, :, :]
        tar_nearWeight = tar_nearWeight[np.newaxis, :, :]
        tar_predictor = tar_predictor[np.newaxis, :, :]


    nstn, ntime = np.shape(stn_data)
    nrow, ncol, nearmax = np.shape(tar_nearIndex)
    estimates = np.nan * np.zeros([nrow, ncol, ntime], dtype=np.float32)

    for r, c in itertools.product(range(nrow), range(ncol)):

        # prepare xdata and sample weight for training and weights of neighboring stations
        sample_nearIndex = tar_nearIndex[r, c, :]
        index_valid = sample_nearIndex >= 0

        if np.sum(index_valid) > 0:
            sample_nearIndex = sample_nearIndex[index_valid]

            sample_weight = tar_nearWeight[r, c, :][index_valid]
            sample_weight = sample_weight / np.sum(sample_weight)

            xdata_near = stn_predictor[sample_nearIndex, :]
            xdata_g = tar_predictor[r, c, :]

            # interpolation for every time step
            for d in range(ntime):

                ydata_near = np.squeeze(stn_data[sample_nearIndex, d])
                if len(np.unique(ydata_near)) == 1:  # e.g., for prcp, all zero
                    ydata_tar = ydata_near[0]
                else:
                    if method == 'linear':
                        ydata_tar = weight_linear_regression(xdata_near, sample_weight, ydata_near, xdata_g)
                    elif method == 'logistic':
                        ydata_tar = weight_logistic_regression(xdata_near, sample_weight, ydata_near, xdata_g)
                    else:
                        sys.exit(f'Unknonwn regression method: {method}')
                estimates[r, c, d] = ydata_tar

    return np.squeeze(estimates)


########################################################################################################################

def main_regression(config, target):
    # target: loo (leave one out station) or grid

    # parse and change configurations
    outpath_parent = config['outpath_parent']
    path_regression = f'{outpath_parent}/regression_outputs'
    os.makedirs(path_regression, exist_ok=True)

    if target == 'grid':
        outfile = f'{path_regression}/Grid_Regression.nc'  # leave one out regression
        config['file_grid_reg'] = outfile
        overwrite_flag = config['overwrite_grid_reg']
        predictor_name_static_target = config['predictor_name_static_grid']
    elif target == 'loo':
        outfile = f'{path_regression}/leave_one_out_Station_Regression.nc'  # leave one out regression
        config['file_loo_reg'] = outfile
        overwrite_flag = config['overwrite_loo_reg']
        predictor_name_static_target = config['predictor_name_static_stn']
    else:
        sys.exit('Unknown target!')

    # in/out information to this function
    file_allstn = config['file_allstn']
    file_stn_nearinfo = config['file_stn_nearinfo']
    file_stn_weight = config['file_stn_weight']
    outfile = outfile # just to make sure all in/out settings are in this section

    target_vars = config['target_vars']
    date_start = config['date_start']
    date_end = config['date_end']
    predictor_name_static_stn = config['predictor_name_static_stn']
    predictor_name_static_target = predictor_name_static_target

    minRange_vars = config['minRange_vars'].copy()
    maxRange_vars = config['maxRange_vars'].copy()
    transform_vars = config['transform_vars']
    transform_settings = config['transform']

    overwrite_flag = overwrite_flag

    if target == 'loo':
        # keyword for near information (default setting in this script)
        near_keyword = 'InStn' # input stations
    else:
        near_keyword = 'Grid'

    print('#' * 50)
    print(f'{target} regression')
    print('#' * 50)
    print('Input file_allstn:', file_allstn)
    print('Input file_stn_nearinfo:', file_stn_nearinfo)
    print('Input file_stn_weight:', file_stn_weight)
    print('Output regression file:', outfile)
    print('Output target:', target)
    print('Target variables:', target_vars)

    if os.path.isfile(outfile):
        print('Note! Output regression file exists')
        if overwrite_flag == True:
            print('overwrite_flag is True. Continue.')
        else:
            print('overwrite_flag is False. Skip regression.')
            return config

    ########################################################################################################################
    # initialize outputs
    with xr.open_dataset(file_allstn) as ds_stn:
        ds_stn = ds_stn.sel(time=slice(date_start, date_end))
        timeaxis = ds_stn.time.values

    with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
        xaxis = ds_nearinfo.x
        yaxis = ds_nearinfo.y

    ds_out = xr.Dataset()
    ds_out.coords['time'] = timeaxis
    if target == 'grid':
        ds_out.coords['x'] = xaxis
        ds_out.coords['y'] = yaxis
    elif target == 'loo':
        ds_out.coords['stn'] = ds_stn.stn.values

    ########################################################################################################################
    # loop variables
    for vn in range(len(target_vars)):

        var_name = target_vars[vn]
        print('Regression for:', var_name)

        # transformed or not
        if len(transform_vars[vn]) > 0:
            var_name_trans = var_name + '_' + transform_vars[vn]
            print(f'Variable {var_name} is transformed using {transform_vars[vn]}')
            print(f'{var_name_trans} instead of {var_name} will be loaded from the station data file {file_allstn}.')

            # adjust max/min limits
            minRange_vars[vn] = data_transformation(minRange_vars[vn], transform_vars[vn], transform_settings[transform_vars[vn]], 'transform')
            maxRange_vars[vn] = data_transformation(maxRange_vars[vn], transform_vars[vn], transform_settings[transform_vars[vn]], 'transform')

        else:
            var_name_trans = ''

        ########################################################################################################################
        # load data for regression

        # station data
        with xr.open_dataset(file_allstn) as ds_stn:
            ds_stn = ds_stn.sel(time=slice(date_start, date_end))

            if len(var_name_trans) > 0:
                stn_value = ds_stn[var_name_trans].values
            else:
                stn_value = ds_stn[var_name].values

            nstn = len(ds_stn.stn)
            predictor_static_stn = np.ones([nstn, len(predictor_name_static_stn) + 1])  # first column used for regression
            for i in range(len(predictor_name_static_stn)):
                predictor_static_stn[:, i + 1] = ds_stn[predictor_name_static_stn[i]].values

        # near information
        with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
            # nearDistance = ds_nearinfo['nearDistance_InStn_' + var_name].values
            vtmp = f'nearIndex_{near_keyword}_{var_name}'
            if vtmp in ds_nearinfo.data_vars:
                nearIndex = ds_nearinfo[vtmp].values
            else:
                if var_name == 'trange':
                    print(f'Use nearIndex_{near_keyword}_tmean for trange')
                    nearIndex = ds_nearinfo[f'nearIndex_{near_keyword}_tmean'].values
                else:
                    sys.exit(f'Cannot find nearIndex_{near_keyword}_{var_name} in {file_stn_nearinfo}')

            if target == 'loo':
                nearIndex = nearIndex[np.newaxis, :, :]

        # predictor information
        if target == 'grid':
            with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
                nrow, ncol, nearmax = np.shape(nearIndex)
                predictor_static_target = np.ones([nrow, ncol, len(predictor_name_static_stn) + 1])  # first column used for regression
                for i in range(len(predictor_name_static_target)):
                    prei = ds_nearinfo[predictor_name_static_target[i]].values
                    predictor_static_target[:, :, i + 1] = prei
        else:
            predictor_static_target = predictor_static_stn.copy()
            predictor_static_target = predictor_static_target[np.newaxis, :, :]

        # weights
        with xr.open_dataset(file_stn_weight) as ds_weight:
            vtmp = f'nearWeight_{near_keyword}_{var_name}'
            if vtmp in ds_weight.data_vars:
                nearWeight = ds_weight[vtmp].values
            else:
                if var_name == 'trange':
                    print(f'Use nearWeight_{near_keyword}_tmean for trange')
                    nearWeight = ds_weight[f'nearWeight_{near_keyword}_tmean'].values
                else:
                    sys.exit(f'Cannot find nearIndex_{near_keyword}_{var_name} in {file_stn_weight}')

            if target == 'loo':
                nearWeight = nearWeight[np.newaxis, :, :]


        ########################################################################################################################
        # produce predictor matrix for regression
        # other static or dynamic predictors can be added in the future

        stn_predictor = predictor_static_stn
        tar_predictor = predictor_static_target
        del predictor_static_stn, predictor_static_target

        ########################################################################################################################
        # get estimates at station points
        estimates = loop_regression_2Dor3D(stn_value, stn_predictor, nearIndex, nearWeight, tar_predictor, 'linear')

        # constrain trange
        estimates = np.squeeze(estimates)
        if np.any(estimates < minRange_vars[vn]):
            print(f'{var_name} estimates have values < {minRange_vars[vn]}. Adjust those to {minRange_vars[vn]}.')
            estimates[estimates < minRange_vars[vn]] = minRange_vars[vn]

        if np.any(estimates > maxRange_vars[vn]):
            print(f'{var_name} estimates have values > {maxRange_vars[vn]}. Adjust those to {maxRange_vars[vn]}.')
            estimates[estimates > maxRange_vars[vn]] = maxRange_vars[vn]

        ########################################################################################################################
        # add to output ds
        if len(var_name_trans) > 0:
            var_name_save = var_name_trans
        else:
            var_name_save = var_name

        if estimates.ndim == 3:
            ds_out[var_name_save] = xr.DataArray(estimates, dims=('y', 'x', 'time'))
        elif estimates.ndim == 2:
            ds_out[var_name_save] = xr.DataArray(estimates, dims=('stn', 'time'))

        ########################################################################################################################
        # if the variable is prcp, do logistic regression too
        if var_name == 'prcp':
            print('Add pop logistic regression')
            if len(var_name_trans) > 0: # this means previous step uses transformed precipitation, while for logistic regression, we use raw precipitation
                stn_value = ds_stn[var_name].values
                print('negative values', np.sum(stn_value<0))
            stn_value[stn_value > 0] = 1
            estimates = loop_regression_2Dor3D(stn_value, stn_predictor, nearIndex, nearWeight, tar_predictor, 'logistic')
            if estimates.ndim == 3:
                ds_out['pop'] = xr.DataArray(estimates, dims=('y', 'x', 'time'))
            elif estimates.ndim == 2:
                ds_out['pop'] = xr.DataArray(estimates, dims=('stn', 'time'))

    # save output file
    encoding = {}
    for var in ds_out.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}

    ds_out.to_netcdf(outfile, encoding=encoding)

    return config
