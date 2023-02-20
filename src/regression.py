
import numpy as np
import pandas as pd
import xarray as xr
# import multiprocessing
from multiprocessing import Pool
from tqdm.contrib import itertools
from data_processing import data_transformation
import sys, os, math

########################################################################################################################
# ludcmp, lubksb, and linearsolver

# Source: https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/ludcmp.py
# Similar to Fortran GMET solution

# "Numerical Recipes" ludcmp() C code on page 46 translated into Python
# I make it as simple as possible, disregard efficiency.
# Here a is a list of list, n is integer (size of the matrix)
# index is a list, and d is also a list of size 1
# Python list index starts from 0.  So matrix index is from 0 to n-1.
#
import math
import time

def ludcmp(a, n, indx, d):
    d[0] = 1.0
    # looking for the largest a in each row and store it in vv as inverse
    # We need a new list same size as indx, for this we use .copy()
    vv = indx.copy()
    for i in range(0, n):
        big = 0.0
        for j in range(0, n):
            temp = math.fabs(a[i][j])
            if (temp > big):
                big = temp
        vv[i] = 1.0 / big
    #
    # run Crout's algorithm
    for j in range(0, n):
        # top half & bottom part are combined
        # but the upper limit l for k sum is different
        big = 0.0
        for i in range(0, n):
            if (i < j):
                l = i
            else:
                l = j
            sum = a[i][j]
            for k in range(0, l):
                sum -= a[i][k] * a[k][j]
            # end for k
            a[i][j] = sum
            # for bottom half, we keep track which row is larger
            if (i >= j):
                dum = vv[i] * math.fabs(sum)
                if (dum >= big):
                    big = dum
                    imax = i
            # end if (i>= ...)
        # end for i
        # pivoting part, swap row j with row imax, a[j] is a whole row
        if (j != imax):
            dum = a[imax]
            a[imax] = a[j]
            a[j] = dum
            d[0] = - d[0]
            vv[imax] = vv[j]
        # end if (j != ...)
        # divide by the beta diagonal value
        indx[j] = imax
        dum = 1.0 / a[j][j]
        for i in range(j + 1, n):
            a[i][j] *= dum
        # end for i
    # end for j


# end of def ludcmp

# We do backward substitution in lubksb() take the row swapped LU decomposed
# a, size n, and swapping indx, and b vector as input.  The output is
# in b after calling.
def lubksb(a, n, indx, b):
    ii = -1
    # forward
    for i in range(0, n):
        ip = indx[i]
        sum = b[ip]
        b[ip] = b[i]
        if (ii != -1):
            for j in range(ii, i):
                sum -= a[i][j] * b[j]
        elif (sum != 0):
            ii = i
        b[i] = sum
    # bote alpha_{ii} is 1 above
    #  backward
    for i in range(n - 1, -1, -1):
        sum = b[i]
        for j in range(i + 1, n):
            sum -= a[i][j] * b[j]
        b[i] = sum / a[i][i]


# end lubksb()


# unfortunately a is destroyed (become swapped LU)
def linearsolver(a, n, b):
    indx = list(range(n))
    d = [1]
    ludcmp(a, n, indx, d)
    x = b.copy()
    lubksb(a, n, indx, x)
    # print("x=", x)
    return x

# if __name__ == "__main__":
#     # using ludcmp and lubksb
#     t1 = time.time()
#     for i in range(10000):
#         a = [[1, 3, 3, -5], [2, -4, 7, -1], [7, 0.5, 3, -6], [9, -2, 3, 8]]
#         n = 4
#         b = [0, 2, 3, -10]
#         x2 = linearsolver(a, n, b)
#     print(x2)
#     t2 = time.time()
#     print(t2-t1)
#
#     # just using numpy
#     from numpy.linalg import inv
#     import numpy as np
#     t1 = time.time()
#     for i in range(10000):
#         a = np.array([[1, 3, 3, -5], [2, -4, 7, -1], [7, 0.5, 3, -6], [9, -2, 3, 8]])
#         b = np.array([0, 2, 3, -10])
#         ainv = inv(a)
#         x2 = np.matmul(ainv, b)
#     print(x2)
#     t2 = time.time()
#     print(t2-t1)


########################################################################################################################
# basic regression utility functions

def least_squares_numpy(x, y, tx):
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


def least_squares_ludcmp(x, y, tx):
    # In fortran version, ludcmp and lubksb are used to calcualte matrix inversion
    # call ludcmp(a, indx, d)
    # call lubksb(a, indx, b)

    # In Python version, numpy is used to calculate matrix inversion
    c = np.matmul(tx, y)
    a = np.matmul(tx, x)
    n = np.shape(a)[0]
    b = linearsolver(list(a), n, list(c))

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
            bn = least_squares_ludcmp(xv, yn, tx)

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
    b = least_squares_ludcmp(nearinfo, datanear, twx_red)
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

def check_predictor_matrix_behavior(nearinfo, weightnear):
    # check to see if the station predictor matrix will be well behaved
    w_pcp_red = np.diag(np.squeeze(weightnear))
    tx_red = np.transpose(nearinfo)
    twx_red = np.matmul(tx_red, w_pcp_red)
    tmp = np.matmul(twx_red, nearinfo)
    vv =  np.max(np.abs(tmp), axis=1)
    if np.any(vv==0):
        flag = False # bad performance
    else:
        flag = True # good performance
    return flag


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
# dynmaic predictor-related functions
def initial_check_dynamic_predictor(dynamic_predictor_name, dynamic_predictor_filelist, target_vars):

    if not isinstance(dynamic_predictor_name, list):
        sys.exit(f'Error! dynamic_predictor_name must be a list.')

    if (len(target_vars) != len(dynamic_predictor_name)) and (len(dynamic_predictor_name) > 0):
        sys.exit(f'Error! len(dynamic_predictor_name)>0 but is different from len(target_vars)!')

    flag = False
    if not os.path.isfile(dynamic_predictor_filelist):
        print('Do not find dynamic_predictor_filelist:', dynamic_predictor_filelist)
    elif len(dynamic_predictor_name) == 0:
        print('dynamic_predictor_name length is 0')
    else:
        with open(dynamic_predictor_filelist, 'r') as f:
            file0 = f.readline().strip()
        if not os.path.isfile(file0):
            print(f'Do not find the first file: {file0} in dynamic_predictor_filelist: {dynamic_predictor_filelist}')
        else:
            # with nc.Dataset(file0) as ds:
            with xr.open_dataset(file0) as ds:
                for i in range(len(dynamic_predictor_name)):
                    tmp = []
                    for v in dynamic_predictor_name[i]:
                        # if v in ds.variables.keys():
                        if v in ds.data_vars:
                            print(f'Include {v} as a dynamic predictor for {target_vars[i]}')
                            tmp.append(v)
                        else:
                            print(f'Cannot find {v} in {file0}. Do not include it as a dynamic predictor for {target_vars[i]}')
                    if len(tmp) > 0:
                        flag = True
                    dynamic_predictor_name[i] = tmp

    if flag == False:
        print('Dynamic predictor is not activated')
    else:
        print(f'Dynamic predictor is activated. Dynamic predicotrs are {dynamic_predictor_name}')
    return flag

def map_filelist_timestep(dynamic_predictor_filelist, timeseries):
    # for every time step, find the corresponding input file

    # get file list
    filelist = []
    with open(dynamic_predictor_filelist, 'r') as f:
        for line in f:
            filelist.append(line.strip())

    # make a dateframe containing all time steps
    df = pd.DataFrame()
    for i in range(len(filelist)):
        with xr.open_dataset(filelist[i]) as ds:
            timei = ds.time.values
            files = [filelist[i]] * len(timei)
            dfi = pd.DataFrame({'intime': timei, 'file': files, 'ind': np.arange(len(timei))})
            df = pd.concat((df, dfi))

    # for each time step, find the closest
    df_mapping = pd.DataFrame()
    for t in timeseries:
        if t >= df['intime'].iloc[0] and t <= df['intime'].iloc[-1]: # don't extrapolate 
            indi = np.argmin(np.abs(df['intime'] - t)) # nearest neighbor
            df_mapping = pd.concat((df_mapping, df.iloc[[indi]]))
        else:
            df_mapping = pd.concat((df_mapping, pd.DataFrame({'intime': [np.nan], 'file': [''], 'ind': [np.nan]})))

    df_mapping['tartime'] = timeseries
    df_mapping.index = np.arange(len(df_mapping))

    return df_mapping

# def read_timestep_input_data(df_mapping, tarindex, varnames):
#     # not finished
#     df = df_mapping.iloc[tarindex]
#     with xr.open_dataset(df['file']) as ds:
#         ds = ds.isel(time=tarindex)
#         ds = ds[varnames]
#         ds = ds.load()
#     return ds


def read_period_input_data(df_mapping, varnames):
    # read and
    # select period
    # select variables

    # tarlon/tarlat: target lat/lon vector
    # default: dynmaic input files must have lat/lon dimensions

    files = np.unique(df_mapping['file'].values)
    files = [f for f in files if len(f)>0]

    intime = df_mapping['intime'].values
    intime = intime[~np.isnan(intime)]

    tartime = df_mapping['tartime'].values

    if len(files) == 0:
        print('Warning! Cannot find any valid dynamic input files')
        ds = xr.Dataset()
    else:
        ds = xr.open_mfdataset(files)
        if np.any([not d in ds.dims for d in ['time', 'lat','lon']]):
            sys.exit(f'Dynamic input files must have lat, lon, and time dimensions!')
        # select ...
        ds = ds[varnames]
        ds = ds.sel(time=slice(tartime[0], tartime[-1]))
        ds = ds.load()
        ds = ds.interp(time=tartime, method='nearest')
    return ds


def regrid_xarray(ds, tarlon, tarlat, target):
    # if target='1D', tarlon and tarlat are vector of station points
    # if target='2D', tarlon and tarlat are vector defining grids

    # regridding: space and time
    if target == '2D':

        ds = ds.transpose('time', 'lat', 'lon')

        if ds.lat.values[0] > ds.lat.values[1]:
            ds = ds.sel(lat=slice(tarlat[-1], tarlat[0]))
        else:
            ds = ds.sel(lat=slice(tarlat[0], tarlat[-1]))
        if ds.lon.values[0] > ds.lon.values[1]:
            ds = ds.sel(lon=slice(tarlon[-1], tarlon[0]))
        else:
            ds = ds.sel(lon=slice(tarlon[0], tarlon[-1]))

        ds_out = ds.interp(lat=tarlat, lon=tarlon, method='linear')

    elif target == '1D':
        tarlat = xr.DataArray(tarlat, dims=('z'))
        tarlon = xr.DataArray(tarlon, dims=('z'))
        ds_out = ds.interp(lat=tarlat, lon=tarlon, method='linear')
        ds_out = ds_out.transpose('time', 'z')

    else:
        sys.exit('Unknown target')

    return ds_out



def flatten_list(lst):
    flat_list = []
    for item in lst:
        if isinstance(item, list):
            flat_list.extend(flatten_list(item))
        else:
            flat_list.append(item)
    return flat_list

########################################################################################################################
# wrapped up regression functions

def loop_regression_2Dor3D(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight, tar_predictor, method, dynamic_predictors={}):
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

    # np.savez_compressed('../docs/data_for_parallel_test.npz', stn_data=stn_data, stn_predictor=stn_predictor,
    #                     tar_nearIndex=tar_nearIndex, tar_nearWeight=tar_nearWeight, tar_predictor=tar_predictor, method=method, dynamic_predictors=dynamic_predictors)

    if len(dynamic_predictors) == 0:
        dynamic_predictors['flag'] = False

    # if tar_nearIndex.ndim == 2:
    #     # make it a 3D array to be consistent
    #     tar_nearIndex = tar_nearIndex[np.newaxis, :, :]
    #     tar_nearWeight = tar_nearWeight[np.newaxis, :, :]
    #     tar_predictor = tar_predictor[np.newaxis, :, :]
    #
    #     if dynamic_predictors['flag'] == True:
    #         # change raw dim: [n_feature, n_time, n_station] to [n_feature, n_time, 1, n_station]
    #         dynamic_predictors['tar_predictor_dynamic'] = dynamic_predictors['tar_predictor_dynamic'][:, :, np.newaxis, :]


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

            xdata_near0 = stn_predictor[sample_nearIndex, :]
            xdata_g0 = tar_predictor[r, c, :]

            # interpolation for every time step
            for d in range(ntime):

                ydata_near = np.squeeze(stn_data[sample_nearIndex, d])
                if len(np.unique(ydata_near)) == 1:  # e.g., for prcp, all zero
                    ydata_tar = ydata_near[0]
                else:

                    # add dynamic predictors if flag is true and predictors are good
                    xdata_near = xdata_near0
                    xdata_g = xdata_g0
                    if dynamic_predictors['flag'] == True:
                        xdata_near_add = dynamic_predictors['stn_predictor_dynamic'][:, d, sample_nearIndex].T
                        xdata_g_add = dynamic_predictors['tar_predictor_dynamic'][:, d, r, c]
                        if np.all(~np.isnan(xdata_near_add)) and np.all(~np.isnan(xdata_g_add)):
                            xdata_near_try = np.hstack((xdata_near, xdata_near_add))
                            xdata_g_try = np.hstack((xdata_g, xdata_g_add))
                            # check if dynamic predictors are good for regression
                            if check_predictor_matrix_behavior(xdata_near_try, sample_weight) == True:
                                xdata_near = xdata_near_try
                                xdata_g = xdata_g_try
                            else:
                                xdata_near_try = np.hstack((xdata_near, xdata_near_add[:, ~dynamic_predictors['predictor_checkflag']]))
                                xdata_g_try = np.hstack((xdata_g, xdata_g_add[~dynamic_predictors['predictor_checkflag']]))
                                if check_predictor_matrix_behavior(xdata_near_try, sample_weight) == True:
                                    xdata_near = xdata_near_try
                                    xdata_g = xdata_g_try

                    # regression
                    if method == 'linear':
                        ydata_tar = weight_linear_regression(xdata_near, sample_weight, ydata_near, xdata_g)
                    elif method == 'logistic':
                        ydata_tar = weight_logistic_regression(xdata_near, sample_weight, ydata_near, xdata_g)
                    else:
                        sys.exit(f'Unknonwn regression method: {method}')
                estimates[r, c, d] = ydata_tar
    return np.squeeze(estimates)


########################################################################################################################
# parallel version of loop regression: independent processes and large memory use if there are many cpus

def init_worker(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight, tar_predictor, method, dynamic_predictors):
    # Using a dictionary is not strictly necessary. You can also
    # use global variables.
    global mppool_ini_dict
    mppool_ini_dict = {}
    mppool_ini_dict['stn_data'] = stn_data
    mppool_ini_dict['stn_predictor'] = stn_predictor
    mppool_ini_dict['tar_nearIndex'] = tar_nearIndex
    mppool_ini_dict['tar_nearWeight'] = tar_nearWeight
    mppool_ini_dict['tar_predictor'] = tar_predictor
    mppool_ini_dict['method'] = method
    mppool_ini_dict['dynamic_predictors'] = dynamic_predictors


def regression_for_blocks(r1, r2, c1, c2):
    stn_data = mppool_ini_dict['stn_data']
    stn_predictor = mppool_ini_dict['stn_predictor']
    tar_nearIndex = mppool_ini_dict['tar_nearIndex']
    tar_nearWeight = mppool_ini_dict['tar_nearWeight']
    tar_predictor = mppool_ini_dict['tar_predictor']
    method = mppool_ini_dict['method']
    dynamic_predictors = mppool_ini_dict['dynamic_predictors']

    nstn, ntime = np.shape(stn_data)
    ydata_tar = np.nan * np.zeros([r2-r1, c2-c1, ntime])

    for r in range(r1, r2):
        for c in range(c1, c2):
            # prepare xdata and sample weight for training and weights of neighboring stations
            sample_nearIndex = tar_nearIndex[r, c, :]
            index_valid = sample_nearIndex >= 0

            if np.sum(index_valid) > 0:
                sample_nearIndex = sample_nearIndex[index_valid]

                sample_weight = tar_nearWeight[r, c, :][index_valid]
                sample_weight = sample_weight / np.sum(sample_weight)

                xdata_near0 = stn_predictor[sample_nearIndex, :]
                xdata_g0 = tar_predictor[r, c, :]

                # interpolation for every time step
                for d in range(ntime):

                    ydata_near = np.squeeze(stn_data[sample_nearIndex, d])
                    if len(np.unique(ydata_near)) == 1:  # e.g., for prcp, all zero
                        ydata_tar[r-r1, c-c1, d] = ydata_near[0]
                    else:

                        # add dynamic predictors if flag is true and predictors are good
                        xdata_near = xdata_near0
                        xdata_g = xdata_g0
                        if dynamic_predictors['flag'] == True:
                            xdata_near_add = dynamic_predictors['stn_predictor_dynamic'][:, d, sample_nearIndex].T
                            xdata_g_add = dynamic_predictors['tar_predictor_dynamic'][:, d, r, c]
                            if np.all(~np.isnan(xdata_near_add)) and np.all(~np.isnan(xdata_g_add)):
                                xdata_near_try = np.hstack((xdata_near, xdata_near_add))
                                xdata_g_try = np.hstack((xdata_g, xdata_g_add))
                                # check if dynamic predictors are good for regression
                                if check_predictor_matrix_behavior(xdata_near_try, sample_weight) == True:
                                    xdata_near = xdata_near_try
                                    xdata_g = xdata_g_try
                                else:
                                    xdata_near_try = np.hstack(
                                        (xdata_near, xdata_near_add[:, ~dynamic_predictors['predictor_checkflag']]))
                                    xdata_g_try = np.hstack((xdata_g, xdata_g_add[~dynamic_predictors['predictor_checkflag']]))
                                    if check_predictor_matrix_behavior(xdata_near_try, sample_weight) == True:
                                        xdata_near = xdata_near_try
                                        xdata_g = xdata_g_try

                        # regression
                        if method == 'linear':
                            ydata_tar[r-r1, c-c1, d] = weight_linear_regression(xdata_near, sample_weight, ydata_near, xdata_g)
                        elif method == 'logistic':
                            ydata_tar[r-r1, c-c1, d] = weight_logistic_regression(xdata_near, sample_weight, ydata_near, xdata_g)
                        else:
                            sys.exit(f'Unknonwn regression method: {method}')

    return ydata_tar

def loop_regression_2Dor3D_multiprocessing(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight, tar_predictor, method, dynamic_predictors={}, num_processes=4):

    if len(dynamic_predictors) == 0:
        dynamic_predictors['flag'] = False

    nstn, ntime = np.shape(stn_data)
    nrow, ncol, nearmax = np.shape(tar_nearIndex)

    # parallel regression
    # items = [(r, r + 1, c, c + 1) for r in range(nrow) for c in range(ncol)] # range(r, r+1), range(c, c+1)
    items = [(r, r + 1, 0, ncol) for r in range(nrow)]

    with Pool(processes=num_processes, initializer=init_worker, initargs=(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight, tar_predictor, method, dynamic_predictors)) as pool:
        result = pool.starmap(regression_for_blocks, items)

    # fill estimates to matrix
    estimates = np.nan * np.zeros([nrow, ncol, ntime], dtype=np.float32)
    for i in range(len(items)):
        indi = items[i]
        estimates[indi[0]:indi[1], indi[2]:indi[3], :] = result[i]

    return np.squeeze(estimates)

########################################################################################################################

def main_regression(config, target):
    # target: loo (leave one out station) or grid
    t1 = time.time()

    # parse and change configurations
    outpath_parent = config['outpath_parent']
    path_regression = f'{outpath_parent}/regression_outputs'
    os.makedirs(path_regression, exist_ok=True)

    datestamp = f"{config['date_start'].replace('-', '')}-{config['date_end'].replace('-', '')}"
    if target == 'grid':
        outfile = f'{path_regression}/Grid_Regression_{datestamp}.nc'  # leave one out regression
        config['file_grid_reg'] = outfile
        overwrite_flag = config['overwrite_grid_reg']
        predictor_name_static_target = config['predictor_name_static_grid']
    elif target == 'loo':
        outfile = f'{path_regression}/Leave_one_out_Regression_{datestamp}.nc'  # leave one out regression
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

    dynamic_predictor_name = config['dynamic_predictor_name']
    dynamic_predictor_filelist = config['dynamic_predictor_filelist']
    dynamic_predictor_check = config['dynamic_predictor_check']

    num_processes = config['num_processes']

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
    print('Number of processes:', num_processes)

    if os.path.isfile(outfile):
        print('Note! Output regression file exists')
        if overwrite_flag == True:
            print('overwrite_flag is True. Continue.')
        else:
            print('overwrite_flag is False. Skip regression.')
            return config

    dynamic_flag = initial_check_dynamic_predictor(dynamic_predictor_name, dynamic_predictor_filelist, target_vars)

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
    # load dynamic factors if dynamic_flag == True
    # this costs more memory than reading for each time step, but reduces interpolation time. this could be a challenge for large domain

    if dynamic_flag == True:

        allvars = flatten_list(dynamic_predictor_name)

        df_mapping = map_filelist_timestep(dynamic_predictor_filelist, timeaxis)
        ds_dynamic = read_period_input_data(df_mapping, allvars)

        ds_dynamic_stn = regrid_xarray(ds_dynamic, ds_stn.lon.values, ds_stn.lat.values, '1D')
        if target == 'grid':
            ds_dynamic_tar = regrid_xarray(ds_dynamic, xaxis, yaxis, '2D')
        elif target == 'loo':
            ds_dynamic_tar = ds_dynamic_stn.copy()


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

        # dynmaic predictors
        predictor_dynamic = {}
        predictor_dynamic['flag'] = dynamic_flag

        if dynamic_flag == True:
            # stn_predictor_dynamic dim: [n_feature, n_time, n_station]
            # tar_predictor_dynamic dim: [n_feature, n_time, n_station] or [n_feature, n_time, n_row, n_col]
            predictor_dynamic['stn_predictor_dynamic'] = np.stack([ds_dynamic_stn[v].values for v in dynamic_predictor_name[vn]], axis=0)
            predictor_dynamic['tar_predictor_dynamic'] = np.stack([ds_dynamic_tar[v].values for v in dynamic_predictor_name[vn]], axis=0)
            predictor_dynamic['predictor_checkflag'] = np.array([v in dynamic_predictor_check for v in dynamic_predictor_name[vn]])

            if target == 'loo':
                # change raw dim: [n_feature, n_time, n_station] to [n_feature, n_time, 1, n_station]
                predictor_dynamic['tar_predictor_dynamic'] = predictor_dynamic['tar_predictor_dynamic'][:, :, np.newaxis, :]

        ########################################################################################################################
        # get estimates at station points
        estimates = loop_regression_2Dor3D_multiprocessing(stn_value, stn_predictor, nearIndex, nearWeight, tar_predictor, 'linear', predictor_dynamic, num_processes)

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
                print('Number of negative values', np.sum(stn_value<0))
            stn_value[stn_value > 0] = 1
            estimates = loop_regression_2Dor3D_multiprocessing(stn_value, stn_predictor, nearIndex, nearWeight, tar_predictor, 'logistic', predictor_dynamic, num_processes)
            if estimates.ndim == 3:
                ds_out['pop'] = xr.DataArray(estimates, dims=('y', 'x', 'time'))
            elif estimates.ndim == 2:
                ds_out['pop'] = xr.DataArray(estimates, dims=('stn', 'time'))

    # save output file
    encoding = {}
    for var in ds_out.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}

    ds_out.to_netcdf(outfile, encoding=encoding)

    t2 = time.time()
    print('Time cost (seconds):', t2-t1)
    print('Successful regression!\n\n')

    return config
