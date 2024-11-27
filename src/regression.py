import sys, os, math, pathlib, time
import numpy as np
import pandas as pd
import xarray as xr
from multiprocessing import Pool
from data_processing import data_transformation, calculate_monthly_cdfs
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from evaluate import evaluate_allpoint
from scipy.interpolate import griddata, RegularGridInterpolator
from sklearn.model_selection import KFold
# from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
# from sklearn.neural_network import MLPRegressor, MLPClassifier
import sklearn
from sklearn import *

########################################################################################################################
# ludcmp, lubksb, and linearsolver

# Source: https://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/ludcmp.py
# Similar to Fortran GMET solution

# "Numerical Recipes" ludcmp() C code on page 46 translated into Python
# I make it as simple as possible, disregard efficiency.
# Here a is a list of list, n is integer (size of the matrix)
# index is a list, and d is also a list of size 1
# Python list index starts from 0.  So matrix index is from 0 to n-1.

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
    p = np.zeros(nstn)  # estimated probability of occurrence
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
            yn = yp - p  # difference between station occurrence (0/1) and estimated probability: Bnew -Bold in Martyn and Slater 2006
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

# ### GMET original linear logistic regression
# def weight_linear_regression(nearinfo, weightnear, datanear, tarinfo):
#     # # nearinfo: predictors from neighboring stations
#     # # [station number, predictor number + 1] array with the first column being ones
#     # nearinfo = np.zeros([nnum, npred+1])
#     #
#     # # weightnear: weight of neighboring stations
#     # # [station number, station number] array with weights located in the diagonal
#     # weightnear = np.zeros([nnum, nnum])
#     # for i in range(nnum):
#     #     weightnear[i, i] = 123
#     #
#     # # tarinfo:  predictors from target stations
#     # # [predictor number + 1] vector with the first value being one
#     # tarinfo = np.zeros(npred+1)
#     #
#     # # datanear: data from neighboring stations. [station number] vector
#     # datanear = np.zeros(nnum)

#     # start regression
#     w_pcp_red = np.diag(np.squeeze(weightnear))
#     tx_red = np.transpose(nearinfo)
#     twx_red = np.matmul(tx_red, w_pcp_red)
#     b = least_squares_ludcmp(nearinfo, datanear, twx_red)
#     datatar = np.dot(tarinfo, b)

#     return datatar

# def weight_logistic_regression(nearinfo, weightnear, datanear, tarinfo):

#     try:
#         w_pcp_red = np.diag(np.squeeze(weightnear))

#         # if len(np.unique(datanear)) == 1:
#         #     poe = datanear[0] # all 0 (no rain) or all 1 (rain everywhere)
#         # else:
#         tx_red = np.transpose(nearinfo)
#         twx_red = np.matmul(tx_red, w_pcp_red)
#         b = logistic_regression(nearinfo, twx_red, datanear)
#         if np.all(b == 0) or np.any(np.isnan(b)):
#             poe = np.dot(weightnear, datanear)
#         else:
#             zb = - np.dot(tarinfo, b)
#             poe = 1 / (1 + np.exp(zb))

#     except:
#         poe = sklearn_weight_logistic_regression(nearinfo, weightnear, datanear, tarinfo)

#     return poe
# ### GMET original linear logistic regression


# ### sklearn linear and logistic regression
# def weight_logistic_regression(nearinfo, weightnear, datanear, tarinfo):

#     try:
#         poe = sklearn_weight_logistic_regression(nearinfo, weightnear, datanear, tarinfo)
#     except:
#         poe = np.nan

#     return poe

# def weight_linear_regression(nearinfo, weightnear, datanear, tarinfo):
#     model = LinearRegression()
#     model = model.fit(nearinfo, datanear, sample_weight=weightnear)
#     datatar = model.predict(tarinfo[np.newaxis, :])
#     return datatar
# ### sklearn linear and logistic regression

### statsmodel linear and logistic regression
import statsmodels.api as sm
def weight_logistic_regression(nearinfo, weightnear, datanear, tarinfo):
    try:
        model = sm.Logit(datanear, nearinfo, weights=weightnear).fit(disp=0)
        poe = model.predict(tarinfo) # Returns probability estimates
    except:
        poe = sklearn_weight_logistic_regression(nearinfo, weightnear, datanear, tarinfo)
    return poe

def weight_linear_regression(nearinfo, weightnear, datanear, tarinfo):
    model = sm.WLS(datanear, nearinfo, weights=weightnear).fit()
    return model.predict(tarinfo[np.newaxis, :])
### statsmodel linear and logistic regression

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
def sklearn_weight_linear_regression(nearinfo, weightnear, datanear, tarinfo):
    model = LinearRegression()
    model = model.fit(nearinfo, datanear, sample_weight=weightnear)
    datatar = model.predict(tarinfo[np.newaxis, :])
    return datatar


def sklearn_weight_logistic_regression(nearinfo, weightnear, datanear, tarinfo):

    model = LogisticRegression(solver='liblinear')
    model = model.fit(nearinfo, datanear, sample_weight=weightnear)
    poe = model.predict_proba(tarinfo[np.newaxis,:])[0, 1]
    return poe


########################################################################################################################
# dynmaic predictor-related functions

def check_dynamic_filepath(dynamic_predictor_filelist):
    # check whether files in dynamic_predictor_filelist use relative or absolute path
    # if the relative path is used, change it to the absolute path for simplicity
    with open(dynamic_predictor_filelist, 'r') as f:
        file0 = f.readlines()

    file1 = []
    path0 = str(pathlib.Path(os.path.abspath(dynamic_predictor_filelist)).parent)
    flag = False
    for f in file0:
        f = f.strip()
        if not os.path.isabs(f):
            file1.append(f'{path0}/{f}')
            flag = True
        else:
            file1.append(f)

    if flag == True:
        print(f'Changing the relative path in {dynamic_predictor_filelist} to the absolute path')
        with open(dynamic_predictor_filelist, 'w') as f:
            for fi in file1:
                _ = f.write(f'{fi}\n')

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

        # change relative path to absolute path if needed
        check_dynamic_filepath(dynamic_predictor_filelist)

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
        print('Dynamic predictors are not activated')
    else:
        print(f'Dynamic predictors are activated. Dynamic predictors are {dynamic_predictor_name}')
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
        # select ...
        ds = ds[varnames]
        ds = ds.sel(time=slice(tartime[0], tartime[-1]))
        ds = ds.load()
        ds = ds.interp(time=tartime, method='linear')
    return ds


def regrid_xarray(ds, tarlon, tarlat, target, method):
    # if target='1D', tarlon and tarlat are vector of station points
    # if target='2D', tarlon and tarlat are vector defining grids

    default_method = 'nearest'
    if len(method) == 0:
        method = default_method

    ds = ds.transpose('time', 'lat', 'lon')

    # regridding: space and time
    if target == '2D':

        if tarlat.ndim == 2:
            # check whether this is regular grids
            if np.all(tarlat[:, 0] == tarlat[:, 1]):
                tarlat = tarlat[:, 0]
            elif np.all(tarlat[0, :] == tarlat[1, :]):
                tarlat = tarlat[0, :]

            if np.all(tarlon[:, 0] == tarlon[:, 1]):
                tarlon = tarlon[:, 0]
            elif np.all(tarlon[0, :] == tarlon[1, :]):
                tarlon = tarlon[0, :]


        if tarlat.ndim == 1:

            if ds.lat.values[0] > ds.lat.values[1]:
                ds = ds.sel(lat=slice(tarlat[-1], tarlat[0]))
            else:
                ds = ds.sel(lat=slice(tarlat[0], tarlat[-1]))
            if ds.lon.values[0] > ds.lon.values[1]:
                ds = ds.sel(lon=slice(tarlon[-1], tarlon[0]))
            else:
                ds = ds.sel(lon=slice(tarlon[0], tarlon[-1]))

            if isinstance(method, str):
                ds_out = ds.interp(lat=tarlat, lon=tarlon, method=method, kwargs={"fill_value": "extrapolate"})
            else:
                ds_out = xr.Dataset()
                for v in ds.data_vars:
                    if v in method:
                        methodv = method[v]
                    else:
                        methodv = default_method
                    ds_out[v] = ds[v].interp(lat=tarlat, lon=tarlon, method=methodv, kwargs={"fill_value": "extrapolate"})

        elif tarlat.ndim == 2:
            ds_out = xr.Dataset()
            ds_out.coords['lat'] = np.arange(tarlat.shape[0])
            ds_out.coords['lon'] = np.arange(tarlat.shape[1])
            ds_out.coords['time'] = ds.time.values
            ds_out['latitude'] = xr.DataArray(tarlat, dims=('lat', 'lon'))
            ds_out['longitude'] = xr.DataArray(tarlon, dims=('lat', 'lon'))

            points = np.column_stack((tarlat.flatten(), tarlon.flatten()))

            for v in ds.data_vars:
                if v in method:
                    methodv = method[v]
                else:
                    methodv = default_method

                # Create the interpolator function
                dv0 = ds[v].values
                dv = np.nan * np.zeros([dv0.shape[0], tarlat.shape[0], tarlon.shape[1]])
                for i in range(dv0.shape[0]):
                    interpolator = RegularGridInterpolator((ds.lat.values, ds.lon.values), dv0[i, :, :], bounds_error=False, method=methodv)
                    zz_new = interpolator(points)
                    zz_new = np.reshape(zz_new, tarlat.shape)
                    dv[i, :, :] = zz_new

                ds_out[v] = xr.DataArray(dv, dims=('time', 'lat', 'lon'))

        ds_out = ds_out.transpose('time', 'lat', 'lon')

    elif target == '1D':

        # simplest method if ds fully contains tarlon/tarlat
        # but for the testcase, some stations are outside the boundary defined by grid centers although they are still within the grids.
        # this method will create NaN
        tarlatz = xr.DataArray(tarlat, dims=('z'))
        tarlonz = xr.DataArray(tarlon, dims=('z'))
        if isinstance(method, str):
            ds_out = ds.interp(lat=tarlatz, lon=tarlonz, method=method)
        else:
            ds_out = xr.Dataset()
            for v in ds.data_vars:
                if v in method:
                    methodv = method[v]
                else:
                    methodv = default_method
                ds_out[v] = ds[v].interp(lat=tarlatz, lon=tarlonz, method=methodv)

        ds_out = ds_out.transpose('time', 'z')

        # if there is any NaN values along the boundaries, fill them using nearest neighbor method
        latstep = np.abs(tarlat[0] - tarlat[1])
        lonstep = np.abs(tarlon[0] - tarlon[1])
        varlist = [v for v in ds_out.data_vars]

        datamatch = np.nan * np.zeros([len(ds.time), len(tarlat), len(ds.data_vars)])
        for i in range(len(varlist)):
            datamatch[:, :, i] = ds_out[varlist[i]].values

        for i in range(len(tarlon)):
            if np.all(np.isnan(datamatch[:, i, 0])):
                londiff = np.abs(ds.lon.values - tarlon[i])
                latdiff = np.abs(ds.lat.values - tarlat[i])
                if np.min(latdiff) < latstep * 2 and np.min(londiff) < lonstep * 2: # * 2 to allow a minimum extrapolation
                    ind1 = np.argmin(latdiff)
                    ind2 = np.argmin(londiff)
                    for j in range(len(varlist)):
                        datamatch[:, i, j] = ds[varlist[j]].values[:, ind1, ind2]

        for i in range(len(varlist)):
            ds_out[varlist[i]].values = datamatch[:, :, i]

        ds_out = xr.Dataset()
        ds_out.coords['time'] = ds['time'].values
        ds_out.coords['z'] = np.arange(len(tarlon))
        for j in range(len(varlist)):
            ds_out[varlist[j]] = xr.DataArray(datamatch[:, :, j], dims=('time', 'z'))

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
# machine learning regression

def train_and_return_test(Xtrain, ytrain, Xtest, method, probflag, settings = {}, weight=[]):
    indexvalid = ~np.isnan( np.sum(Xtrain, axis=1) + ytrain)

    if Xtest.ndim == 1:
        Xtest = Xtest[np.newaxis, :]

    indexnan_test = np.isnan(np.sum(Xtest, axis=1))
    Xtest = Xtest.copy()
    Xtest[indexnan_test, :] = 0

    Xtrain = Xtrain[indexvalid, :]
    ytrain = ytrain[indexvalid]

    ldict = {'settings': settings, 'method': method}
    exec(f"model = sklearn.{method}(**settings)", globals(), ldict)
    # exec(f"model = {method}(**settings)", globals(), ldict)
    model = ldict['model']

    if len(weight) == 0:
        model.fit(Xtrain, ytrain)
    else:
        model.fit(Xtrain, ytrain, weight)

    if probflag == False:
        ytest = model.predict(Xtest)
    else:
        try:
            ytest = model.predict_proba(Xtest)[:, 1]
        except:
            ytest = model.predict(Xtest)

    ytest[indexnan_test] = np.nan
    return ytest

def ML_regression_crossvalidation(stn_data, stn_predictor, ml_model, probflag, ml_settings={}, dynamic_predictors={}, n_splits=10):
    t1 = time.time()

    if len(dynamic_predictors) == 0:
        dynamic_predictors['flag'] = False

    # remove missing stations
    nstn0, ntime  = np.shape(stn_data)
    estimates0    = np.nan * np.zeros([nstn0, ntime], dtype=np.float32)
    nannum        = np.sum(np.isnan(stn_data), axis=1)
    indexmissing  = nannum == stn_data.shape[1]
    stn_data      = stn_data[~indexmissing, :]
    stn_predictor = stn_predictor[~indexmissing, :]

    # cross-validation index
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    kf.get_n_splits(stn_data)

    nstn, ntime = np.shape(stn_data)
    estimates = np.nan * np.zeros([nstn, ntime], dtype=np.float32)

    for i, (train_index, test_index) in enumerate(kf.split(stn_predictor)):
        for t in range(ntime):
            stn_predictor_t = stn_predictor.copy()
            if dynamic_predictors['flag'] == True:
                xdata_add = dynamic_predictors['stn_predictor_dynamic'][:, t, :].T
                stn_predictor_t = np.hstack((stn_predictor_t, xdata_add[~indexmissing, :]))
            ytest = train_and_return_test(stn_predictor_t[train_index, :], stn_data[train_index, t],
                                          stn_predictor_t[test_index, :], ml_model, probflag, ml_settings)
            estimates[test_index, t] = ytest

    estimates0[~indexmissing, :] = estimates

    t2 = time.time()
    print('Regression time cost (s):', t2-t1)
    return np.squeeze(estimates0)

def init_worker_sklearn(stn_data, stn_predictor, tar_predictor, ml_model, probflag, ml_settings, dynamic_predictors, train_index, test_index, indexmissing, randseed=123456789):
    # Using a dictionary is not strictly necessary. You can also
    # use global variables.
    global mppool_ini_dict_sk
    mppool_ini_dict_sk = {}
    mppool_ini_dict_sk['stn_data'] = stn_data
    mppool_ini_dict_sk['stn_predictor'] = stn_predictor
    mppool_ini_dict_sk['tar_predictor'] = tar_predictor
    mppool_ini_dict_sk['ml_model'] = ml_model
    mppool_ini_dict_sk['ml_settings'] = ml_settings
    mppool_ini_dict_sk['probflag'] = probflag
    mppool_ini_dict_sk['train_index'] = train_index
    mppool_ini_dict_sk['test_index'] = test_index
    mppool_ini_dict_sk['dynamic_predictors'] = dynamic_predictors
    mppool_ini_dict_sk['indexmissing'] = indexmissing
    mppool_ini_dict_sk['randseed'] = randseed

def train_and_return_test_cv_timestep(t):
    stn_data = mppool_ini_dict_sk['stn_data']
    stn_predictor = mppool_ini_dict_sk['stn_predictor']
    ml_model = mppool_ini_dict_sk['ml_model']
    ml_settings = mppool_ini_dict_sk['ml_settings']
    probflag = mppool_ini_dict_sk['probflag']
    train_indexall = mppool_ini_dict_sk['train_index']
    test_indexall = mppool_ini_dict_sk['test_index']
    dynamic_predictors = mppool_ini_dict_sk['dynamic_predictors']
    indexmissing = mppool_ini_dict_sk['indexmissing']
    randseed = mppool_ini_dict_sk['randseed']

    np.random.seed(randseed + t)

    nstn, ntime = np.shape(stn_data)
    estimates = np.nan * np.zeros(nstn, dtype=np.float32)

    for train_index, test_index in zip(train_indexall, test_indexall):
        stn_predictor_t = stn_predictor.copy()
        if dynamic_predictors['flag'] == True:
            xdata_add = dynamic_predictors['stn_predictor_dynamic'][:, t, :].T
            stn_predictor_t = np.hstack((stn_predictor_t, xdata_add[~indexmissing, :]))
        ytest = train_and_return_test(stn_predictor_t[train_index, :], stn_data[train_index, t], stn_predictor_t[test_index, :], ml_model, probflag, ml_settings)
        estimates[test_index] = ytest

    return estimates

def ML_regression_crossvalidation_multiprocessing(stn_data, stn_predictor, ml_model, probflag, ml_settings={}, dynamic_predictors={}, n_splits=10, num_processes=1, master_seed=123456789):
    t1 = time.time()

    if len(dynamic_predictors) == 0:
        dynamic_predictors['flag'] = False

    # remove missing stations
    nstn0, ntime = np.shape(stn_data)
    estimates0 = np.nan * np.zeros([nstn0, ntime], dtype=np.float32)

    nannum = np.sum(np.isnan(stn_data), axis=1)
    indexmissing = nannum == stn_data.shape[1]
    stn_data = stn_data[~indexmissing, :]
    stn_predictor = stn_predictor[~indexmissing, :]

    # cross-validation index
    kf = KFold(n_splits=n_splits, shuffle=True)
    kf.get_n_splits(stn_data)
    train_index = []
    test_index = []
    for i, (ind1, ind2) in enumerate(kf.split(stn_predictor)):
        train_index.append(ind1)
        test_index.append(ind2)

    # parallel regression
    items = [(t,) for t in range(ntime)]
    with Pool(processes=num_processes, initializer=init_worker_sklearn, initargs=(stn_data, stn_predictor, [], ml_model, probflag, ml_settings, dynamic_predictors, train_index, test_index, indexmissing, master_seed)) as pool:
        result = pool.starmap(train_and_return_test_cv_timestep, items)

    # fill estimates to matrix
    nstn, ntime = np.shape(stn_data)
    estimates = np.nan * np.zeros([nstn, ntime], dtype=np.float32)
    for i in range(len(items)):
        indi = items[i]
        estimates[:, indi[0]] = result[i]

    estimates0[~indexmissing, :] = estimates

    t2 = time.time()
    print('Regression time cost (sec):', t2-t1)
    return np.squeeze(estimates0)

def ML_regression_grid(stn_data, stn_predictor, tar_predictor, ml_model, probflag, ml_settings={}, dynamic_predictors={}):
    t1 = time.time()

    if len(dynamic_predictors) == 0:
        dynamic_predictors['flag'] = False

    # remove missing stations
    nstn, ntime = np.shape(stn_data)
    nrow, ncol, npred = np.shape(tar_predictor)
    estimates = np.nan * np.zeros([nrow, ncol, ntime], dtype=np.float32)

    tar_predictor_resh = np.reshape(tar_predictor, [nrow*ncol, npred])

    for t in range(ntime):
        stn_predictor_t = stn_predictor.copy()
        tar_predictor_t = tar_predictor_resh.copy()

        if dynamic_predictors['flag'] == True:
            xdata_add1 = dynamic_predictors['stn_predictor_dynamic'][:, t, :].T
            stn_predictor_t = np.hstack((stn_predictor_t, xdata_add1))

            xdata_add2 = dynamic_predictors['tar_predictor_dynamic'][:, t, :, :]
            ndyn = xdata_add2.shape[0]
            xdata_add2 = np.reshape(xdata_add2, [ndyn, nrow*ncol]).T
            tar_predictor_t = np.hstack((tar_predictor_t, xdata_add2))

        ytest = train_and_return_test(stn_predictor_t, stn_data[:, t], tar_predictor_t, ml_model, probflag, ml_settings)
        ytest = np.reshape(ytest, [nrow, ncol])
        estimates[:, :, t] = ytest

    t2 = time.time()
    print('Regression time cost (sec):', t2-t1)
    return np.squeeze(estimates)

def train_and_return_test_grid_timestep(t):
    stn_data = mppool_ini_dict_sk['stn_data']
    stn_predictor = mppool_ini_dict_sk['stn_predictor']
    ml_model = mppool_ini_dict_sk['ml_model']
    ml_settings = mppool_ini_dict_sk['ml_settings']
    probflag = mppool_ini_dict_sk['probflag']
    dynamic_predictors = mppool_ini_dict_sk['dynamic_predictors']
    tar_predictor = mppool_ini_dict_sk['tar_predictor']
    randseed = mppool_ini_dict_sk['randseed']

    np.random.seed(randseed + t)

    nstn, ntime = np.shape(stn_data)
    nrow, ncol, npred = np.shape(tar_predictor)
    tar_predictor_resh = np.reshape(tar_predictor, [nrow * ncol, npred])

    stn_predictor_t = stn_predictor.copy()
    tar_predictor_t = tar_predictor_resh.copy()

    if dynamic_predictors['flag'] == True:
        xdata_add1 = dynamic_predictors['stn_predictor_dynamic'][:, t, :].T
        stn_predictor_t = np.hstack((stn_predictor_t, xdata_add1))

        xdata_add2 = dynamic_predictors['tar_predictor_dynamic'][:, t, :, :]
        ndyn = xdata_add2.shape[0]
        xdata_add2 = np.reshape(xdata_add2, [ndyn, nrow*ncol]).T
        tar_predictor_t = np.hstack((tar_predictor_t, xdata_add2))

    ytest = train_and_return_test(stn_predictor_t, stn_data[:, t], tar_predictor_t, ml_model, probflag, ml_settings)
    ytest = np.reshape(ytest, [nrow, ncol])

    return ytest

def ML_regression_grid_multiprocessing(stn_data, stn_predictor, tar_predictor, ml_model, probflag, ml_settings={}, dynamic_predictors={}, num_processes=1, master_seed=123456789):
    t1 = time.time()

    if len(dynamic_predictors) == 0:
        dynamic_predictors['flag'] = False

    # remove missing stations
    nstn, ntime = np.shape(stn_data)
    nrow, ncol, npred = np.shape(tar_predictor)
    estimates = np.nan * np.zeros([nrow, ncol, ntime], dtype=np.float32)

    # parallel regression
    items = [(t,) for t in range(ntime)]
    with Pool(processes=num_processes, initializer=init_worker_sklearn, initargs=(stn_data, stn_predictor, tar_predictor, ml_model, probflag, ml_settings, dynamic_predictors, [], [], [], master_seed)) as pool:
        result = pool.starmap(train_and_return_test_grid_timestep, items)

    # fill estimates to matrix
    for i in range(len(items)):
        indi = items[i]
        estimates[:, :, indi[0]] = result[i]

    t2 = time.time()
    print('Regression time cost (sec):', t2-t1)
    return np.squeeze(estimates)

########################################################################################################################
# parallel version of loop regression: independent processes and large memory use if there are many cpus

def init_worker(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight, tar_predictor, method, probflag, settings,
                dynamic_predictors, importmodules):
    # Using a dictionary is not strictly necessary. You can also
    # use global variables.
    global mppool_ini_dict
    mppool_ini_dict                       = {}
    mppool_ini_dict['stn_data']           = stn_data
    mppool_ini_dict['stn_predictor']      = stn_predictor
    mppool_ini_dict['tar_nearIndex']      = tar_nearIndex
    mppool_ini_dict['tar_nearWeight']     = tar_nearWeight
    mppool_ini_dict['tar_predictor']      = tar_predictor
    mppool_ini_dict['method']             = method
    mppool_ini_dict['probflag']           = probflag
    mppool_ini_dict['settings']           = settings
    mppool_ini_dict['dynamic_predictors'] = dynamic_predictors

    for im in importmodules:
        if '.' in im:
            im2 = im.split('.')[-1]
            im1 = im.replace('.'+im2, '')
            exec(f"from {im1} import {im2}", globals())
        else:
            exec(f"import {im}", globals())

    # global Ridge
    # from sklearn.linear_model import Ridge


def regression_for_blocks(r1, r2, c1, c2):
    stn_data           = mppool_ini_dict['stn_data']
    stn_predictor      = mppool_ini_dict['stn_predictor']
    tar_nearIndex      = mppool_ini_dict['tar_nearIndex']
    tar_nearWeight     = mppool_ini_dict['tar_nearWeight']
    tar_predictor      = mppool_ini_dict['tar_predictor']
    method             = mppool_ini_dict['method']
    probflag           = mppool_ini_dict['probflag']
    settings           = mppool_ini_dict['settings']
    dynamic_predictors = mppool_ini_dict['dynamic_predictors']

    nstn, ntime = np.shape(stn_data)
    ydata_tar   = np.nan * np.zeros([r2-r1, c2-c1, ntime])

    for r in range(r1, r2):
        for c in range(c1, c2):
            # prepare xdata and sample weight for training and weights of neighboring stations
            sample_nearIndex = tar_nearIndex[r, c, :]
            index_valid = sample_nearIndex >= 0

            if np.sum(index_valid) > 0:
                sample_nearIndex = sample_nearIndex[index_valid]

                sample_weight = tar_nearWeight[r, c, :][index_valid]

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

                                # whether use a dynamic predictor
                                diff = np.max(xdata_near_add, axis=0) - np.min(xdata_near_add, axis=0)
                                tolerance = 1e-10

                                xdata_near_add = xdata_near_add[:, diff > tolerance]
                                xdata_g_add = xdata_g_add[diff > tolerance]

                                if xdata_near_add.size > 0:
                                    xdata_near_try = np.hstack((xdata_near, xdata_near_add))
                                    xdata_g_try = np.hstack((xdata_g, xdata_g_add))

                                    # The below codes using check_predictor_matrix_behavior are not necessary
                                    xdata_near = xdata_near_try
                                    xdata_g = xdata_g_try

                                    # # check if dynamic predictors are good for regression. not necessary after unique value check
                                    # if check_predictor_matrix_behavior(xdata_near_try, sample_weight) == True:
                                    #     xdata_near = xdata_near_try
                                    #     xdata_g = xdata_g_try
                                    # else:
                                    #     xdata_near_try = np.hstack((xdata_near, xdata_near_add[:, ~dynamic_predictors['predictor_checkflag']]))
                                    #     xdata_g_try = np.hstack((xdata_g, xdata_g_add[~dynamic_predictors['predictor_checkflag']]))
                                    #     if check_predictor_matrix_behavior(xdata_near_try, sample_weight) == True:
                                    #         xdata_near = xdata_near_try
                                    #         xdata_g = xdata_g_try

                        # regression
                        if method == 'Linear':
                            ydata_tar[r-r1, c-c1, d] = weight_linear_regression(xdata_near, sample_weight, ydata_near, xdata_g)
                        elif method == 'Logistic':
                            ydata_tar[r-r1, c-c1, d] = weight_logistic_regression(xdata_near, sample_weight, ydata_near, xdata_g)
                        else:
                            ydata_tar[r-r1, c-c1, d] = train_and_return_test(xdata_near, ydata_near, xdata_g, method, probflag,
                                                                             settings, sample_weight)
                        # else:
                        #     sys.exit(f'Unknonwn regression method: {method}')

    return ydata_tar

def loop_regression_2Dor3D_multiprocessing(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight, tar_predictor, method, probflag,
                                           settings, dynamic_predictors={}, num_processes=4, importmodules=[]):
#def loop_regression_2Dor3D_multiprocessing(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight, tar_predictor, method, probflag,
#                                           settings, dynamic_predictors={}, num_processes=4):
    t1 = time.time()

    if len(dynamic_predictors) == 0:
        dynamic_predictors['flag'] = False

    nstn, ntime = np.shape(stn_data)
    nrow, ncol, nearmax = np.shape(tar_nearIndex)

    # parallel regression
    # items = [(r, r + 1, c, c + 1) for r in range(nrow) for c in range(ncol)] # range(r, r+1), range(c, c+1)
    if nrow > ncol:
        items = [(r, r + 1, 0, ncol) for r in range(nrow)]
    else:
        items = [(0, nrow, c, c + 1) for c in range(ncol)]

    with Pool(processes=num_processes, initializer=init_worker, initargs=(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight,
                                                                          tar_predictor, method, probflag, settings, dynamic_predictors,
                                                                          importmodules)) as pool:
    #with Pool(processes=num_processes, initializer=init_worker, initargs=(stn_data, stn_predictor, tar_nearIndex, tar_nearWeight,
    #                                                                      tar_predictor, method, probflag, settings,
    #                                                                      dynamic_predictors)) as pool:
        result = pool.starmap(regression_for_blocks, items)

    # fill estimates to matrix
    estimates = np.nan * np.zeros([nrow, ncol, ntime], dtype=np.float32)
    for i in range(len(items)):
        indi = items[i]
        estimates[indi[0]:indi[1], indi[2]:indi[3], :] = result[i]

    t2 = time.time()
    print('Regression time cost (sec):', t2-t1)
    return np.squeeze(estimates)

########################################################################################################################

def main_regression(config, target):
    # target: loo (leave one out station) or grid
    t1 = time.time()

    # parse and change configurations
    case_name = config['case_name']

    outpath_parent  = config['outpath_parent']
    path_regression = f'{outpath_parent}/regression/'
    os.makedirs(path_regression, exist_ok=True)

    datestamp = f"{config['date_start'].replace('-', '')}-{config['date_end'].replace('-', '')}"

    if 'append_date_to_output_filename' in config:
        append_date_to_output_filename = config['append_date_to_output_filename']
    else:
        append_date_to_output_filename = False

    if target == 'grid':
        if append_date_to_output_filename == True:
            outfile = f'{path_regression}/{case_name}_grid_regression_{datestamp}.nc'  # regression without cross-validation
        else:
            outfile = f'{path_regression}/{case_name}_grid_regression.nc'
        config['file_grid_reg'] = outfile
        if 'overwrite_grid_reg' in config:
            overwrite_flag = config['overwrite_grid_reg']
        else:
            overwrite_flag = False
        predictor_name_static_target = config['predictor_name_static_grid']

    elif target == 'cval':
        if append_date_to_output_filename == True:
            outfile = f'{path_regression}/{case_name}_stn_CV_regression_{datestamp}.nc'  # leave one out regression
        else:
            outfile = f'{path_regression}/{case_name}_stn_CV_regression.nc'
        config['file_cval_reg'] = outfile

        if 'overwrite_stn_cv_reg' in config:
            overwrite_flag = config['overwrite_stn_cv_reg']
        else:
            overwrite_flag = False
        predictor_name_static_target = config['predictor_name_static_stn']
    else:
        sys.exit('Unknown target!')

    # in/out information to this function

    file_allstn       = config['file_allstn']
    file_stn_nearinfo = config['file_stn_nearinfo']
    file_stn_weight   = config['file_stn_weight']
    outfile           = outfile # just to make sure all in/out settings are in this section

    target_vars       = config['target_vars']

    if 'target_vars_WithProbability' in config:
        target_vars_WithProbability = config['target_vars_WithProbability']
    else:
        target_vars_WithProbability = []

    if 'probability_thresholds' in config:
        probability_thresholds      = config['probability_thresholds']
    else:
        probability_thresholds = [0] * len(target_vars_WithProbability)

    date_start        = config['date_start']
    date_end          = config['date_end']

    predictor_name_static_stn    = config['predictor_name_static_stn']
    predictor_name_static_target = predictor_name_static_target

    if 'minRange_vars' in config:
        minRange_vars = config['minRange_vars']
        if not isinstance(minRange_vars, list):
            minRange_vars = [minRange_vars] * len(target_vars)
        minRange_vars = minRange_vars.copy()
    else:
        minRange_vars = [-np.inf] * len(target_vars)

    if 'maxRange_vars' in config:
        maxRange_vars = config['maxRange_vars']
        if not isinstance(maxRange_vars, list):
            maxRange_vars = [maxRange_vars] * len(target_vars)
        maxRange_vars = maxRange_vars.copy()
    else:
        maxRange_vars = [np.inf] * len(target_vars)

    if 'transform_vars' in config:
        transform_vars = config['transform_vars']
        if not isinstance(transform_vars, list):
            transform_vars = [transform_vars] * len(target_vars)
    else:
        transform_vars = [''] * len(target_vars)

    if 'transform' in config:
        transform_settings = config['transform']
    else:
        transform_settings = {}

    if 'dynamic_predictor_filelist' in config:
        dynamic_predictor_filelist  = config['dynamic_predictor_filelist']
        dynamic_predictor_name      = config['dynamic_predictor_name']
        dynamic_predictor_operation = config['dynamic_predictor_operation']
    else:
        dynamic_predictor_filelist  = ''
        dynamic_predictor_name      = []
        dynamic_predictor_operation = []

    num_processes = config['num_processes']
    if 'master_seed' in config:
        master_seed = config['master_seed']
    else:
        master_seed = -1
    overwrite_flag = overwrite_flag

    gridcore_classification = config['gridcore_classification']
    gridcore_continuous     = config['gridcore_continuous']

    n_splits = config['n_splits']

    ensemble_flag = config['ensemble_flag']
    backtransform = not ensemble_flag       # for example, if ensemble_flag=False, no need to create ensemble outputs,
                                            #   the regression outputs should be backtransformed in this step

    if 'sklearn' in config:
        sklearn_config = config['sklearn']
    else:
        sklearn_config = {}

    if target == 'cval':
        # keyword for near information (default setting in this script)
        near_keyword = 'InStn'   # input stations
    else:
        near_keyword = 'Grid'

    stn_lat_name = config['stn_lat_name']
    stn_lon_name = config['stn_lon_name']

    grid_lat_name = config['grid_lat_name']
    grid_lon_name = config['grid_lon_name']

    dynamic_grid_lat_name = config['dynamic_grid_lat_name']
    dynamic_grid_lon_name = config['dynamic_grid_lon_name']

    print('#' * 50)
    if(target == 'cval'):
        print('step 1:  station point cross-validated regression to estimate predictive uncertainty')
    elif(target == 'grid'):
        print('step 2:  regression for all grid points')
    else:
        print('regression target not recognized:', target)
    print('#' * 50)
    print('Input file_allstn:      ', file_allstn)
    print('Input file_stn_nearinfo:', file_stn_nearinfo)
    print('Input file_stn_weight:  ', file_stn_weight)
    print('Output regression file: ', outfile)
    print('Output target:          ', target)
    print('Target variables:       ', target_vars)
    print('Number of processes:    ', num_processes)

    if os.path.isfile(outfile):
        print('Note! Output regression file exists')
        if overwrite_flag == True:
            print('overwrite_flag is True. Continue.')
        else:
            print('overwrite_flag is False. Skip regression.\n')
            return config

    dynamic_flag = initial_check_dynamic_predictor(dynamic_predictor_name, dynamic_predictor_filelist, target_vars)

    if master_seed < 0:
        master_seed = np.random.randint(1e9)
    np.random.seed(master_seed)


    ########################################################################################################################
    # load sklearn methods and settings
    gnew = []
    importmodules = []
    for g in [gridcore_classification, gridcore_continuous]:
        if g.startswith('LWR:'):
            g = g.replace('LWR:', '')

        # load sklearn modules
        if '.' in g:
            m1 = g.split('.')[0]
            m2 = g.split('.')[1]
            # exec(f'global {m2}')
            # exec(f"from sklearn.{m1} import {m2}")
            g = m2
            importmodules.append(f'sklearn.{m1}.{m2}')

        gnew.append(g)

        # sklearn settings
        if not g in sklearn_config:
            sklearn_config[g] = {}

    gridcore_classification_short, gridcore_continuous_short = gnew

    ########################################################################################################################
    # initialize outputs
    with xr.open_dataset(file_allstn) as ds_stn:
        ds_stn = ds_stn.sel(time=slice(date_start, date_end))
        timeaxis = ds_stn.time.values

    ds_out = xr.Dataset()
    ds_out.coords['time'] = timeaxis
    if target == 'grid':
        with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
            xaxis = ds_nearinfo[grid_lon_name].isel(y=0).values
            yaxis = ds_nearinfo[grid_lat_name].isel(x=0).values
            ds_out.coords['x'] = ds_nearinfo.coords['x']
            ds_out.coords['y'] = ds_nearinfo.coords['y']
            ds_out[grid_lat_name] = ds_nearinfo[grid_lat_name]
            ds_out[grid_lon_name] = ds_nearinfo[grid_lon_name]

    elif target == 'cval':
        ds_out.coords['stn'] = ds_stn.stn.values

    else:
        sys.exit(f'Unknown target: {target}')

    ########################################################################################################################
    # load dynamic factors if dynamic_flag == True
    # this costs more memory than reading for each time step, but reduces interpolation time. this could be a challenge for large domain

    if dynamic_flag == True:

        allvars = flatten_list(dynamic_predictor_name)

        df_mapping = map_filelist_timestep(dynamic_predictor_filelist, timeaxis)
        ds_dynamic = read_period_input_data(df_mapping, allvars)
        ds_dynamic = ds_dynamic.rename({dynamic_grid_lat_name:'lat', dynamic_grid_lon_name:'lon'})

        # transformation dynamic variables if necessary
        dyn_operation_trans = {}
        dyn_operation_interp = {}
        for op in dynamic_predictor_operation:
            info = op.split(':')
            for i in range(1, len(info)):
                if info[i].split('=')[0] == 'transform':
                    dyn_operation_trans[info[0]] = info[i].split('=')[1]
                elif info[i].split('=')[0] == 'interp':
                    dyn_operation_interp[info[0]] = info[i].split('=')[1]

        for v in ds_dynamic.data_vars:
            if v in dyn_operation_trans:
                print('Transform dynamic predictor:', v)
                ds_dynamic[v].values = data_transformation(ds_dynamic[v].values, dyn_operation_trans[v],
                                                           transform_settings[dyn_operation_trans[v]], 'transform')

        ds_dynamic_stn = regrid_xarray(ds_dynamic, ds_stn[stn_lon_name].values, ds_stn[stn_lat_name].values, '1D', method=dyn_operation_interp)
        if target == 'grid':
            # ds_dynamic_tar2 = regrid_xarray(ds_dynamic, xaxis, yaxis, '2D', method=dyn_operation_interp)
            ds_dynamic_tar = regrid_xarray(ds_dynamic, ds_nearinfo[grid_lon_name].values, ds_nearinfo[grid_lat_name].values, '2D',
                                           method=dyn_operation_interp)

        elif target == 'cval':
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
            if minRange_vars[vn] != -np.inf:
                minRange_vars[vn] = data_transformation(minRange_vars[vn], transform_vars[vn], transform_settings[transform_vars[vn]], 'transform')
            if maxRange_vars[vn] != np.inf:
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
                sys.exit(f'Cannot find nearIndex_{near_keyword}_{var_name} in {file_stn_nearinfo}')

            if target == 'cval':
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
                sys.exit(f'Cannot find nearIndex_{near_keyword}_{var_name} in {file_stn_weight}')

            if target == 'cval':
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

            if target == 'cval':
                # change raw dim: [n_feature, n_time, n_station] to [n_feature, n_time, 1, n_station]
                predictor_dynamic['tar_predictor_dynamic'] = predictor_dynamic['tar_predictor_dynamic'][:, :, np.newaxis, :]

        ########################################################################################################################
        # get estimates at station points
        probflag = False  # for continuous variables
        if gridcore_continuous.startswith('LWR:'):
            estimates = loop_regression_2Dor3D_multiprocessing(stn_value, stn_predictor, nearIndex, nearWeight, tar_predictor,
                                                               gridcore_continuous[4:], probflag, sklearn_config[gridcore_continuous_short],
                                                               predictor_dynamic, num_processes)
        else:
            if target == 'cval':
                estimates = ML_regression_crossvalidation_multiprocessing(stn_value, stn_predictor, gridcore_continuous, probflag,
                                                                          sklearn_config[gridcore_continuous_short], predictor_dynamic,
                                                                          n_splits, num_processes, master_seed)
            else:
                estimates = ML_regression_grid_multiprocessing(stn_value, stn_predictor, tar_predictor, gridcore_continuous, probflag,
                                                               sklearn_config[gridcore_continuous_short], predictor_dynamic, num_processes,
                                                               master_seed)
        # constrain variables
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
            if backtransform == False:
                var_name_save = var_name_trans
            else:
                var_name_save = var_name
                if 'ecdf' in var_name_trans:
                    cdfs = calculate_monthly_cdfs(xr.open_dataset(file_allstn),var_name,transform_settings[transform_vars[vn]])
                    estimates = data_transformation(estimates,transform_vars[vn], transform_settings[transform_vars[vn]], 'back_transform',times=ds_out['time'].values, cdfs=cdfs)
                else:   
                    estimates = data_transformation(estimates,transform_vars[vn], transform_settings[transform_vars[vn]], 'back_transform')
        else:
            var_name_save = var_name

        if estimates.ndim == 3:
            ds_out[var_name_save] = xr.DataArray(estimates, dims=('y', 'x', 'time'))
        elif estimates.ndim == 2:
            ds_out[var_name_save] = xr.DataArray(estimates, dims=('stn', 'time'))

            # evaluation
            dtmp1 = ds_stn[var_name].values
            if (len(var_name_trans) > 0) and (backtransform == False):
                if 'ecdf' in var_name_trans:
                    cdfs = calculate_monthly_cdfs(xr.open_dataset(file_allstn),var_name,transform_settings[transform_vars[vn]]) 
                    dtmp2 = data_transformation(estimates,transform_vars[vn], transform_settings[transform_vars[vn]], 'back_transform',times=ds_out['time'].values, cdfs=cdfs)
                else:
                    dtmp2 = data_transformation(estimates,transform_vars[vn], transform_settings[transform_vars[vn]], 'back_transform')
            else:
                dtmp2 = estimates

            metvalue, metname = evaluate_allpoint(dtmp1, dtmp2, np.nan)
            ds_out.coords['met'] = metname
            ds_out[var_name_save + '_metric'] = xr.DataArray(metvalue, dims=('stn', 'met'))
            del dtmp1, dtmp2

        ########################################################################################################################
        # if the variable has occurrence features, do logistic regression too
        if var_name in target_vars_WithProbability:
            print(f'Add prob. of event occurrence (POE) for {var_name}: it is in target_vars_WithProbability {target_vars_WithProbability}')
            stn_value = ds_stn[var_name].values.copy()
            var_threshold = probability_thresholds[target_vars_WithProbability.index(var_name)]
            # print(f'Number of values <= threshold {var_threshold}', np.sum(stn_value<=var_threshold))
            stn_value[stn_value <= var_threshold] = 0
            stn_value[stn_value > var_threshold] = 1

            probflag = True
            if gridcore_classification.startswith('LWR:'):

                estimates = loop_regression_2Dor3D_multiprocessing(stn_value, stn_predictor, nearIndex, nearWeight, tar_predictor,
                                                                   gridcore_classification[4:], probflag,
                                                                   sklearn_config[gridcore_continuous_short], predictor_dynamic,
                                                                   num_processes, importmodules)
            else:
                if target == 'cval':
                    estimates = ML_regression_crossvalidation_multiprocessing(stn_value, stn_predictor, gridcore_classification, probflag,
                                                                              sklearn_config[gridcore_classification_short], predictor_dynamic,
                                                                              n_splits, num_processes, master_seed)
                else:
                    estimates = ML_regression_grid_multiprocessing(stn_value, stn_predictor, tar_predictor, gridcore_classification, probflag,
                                                                   sklearn_config[gridcore_classification_short], predictor_dynamic,
                                                                   num_processes, master_seed)

            var_poe = var_name + '_poe'
            if estimates.ndim == 3:
                ds_out[var_poe] = xr.DataArray(estimates, dims=('y', 'x', 'time'))
            elif estimates.ndim == 2:
                ds_out[var_poe] = xr.DataArray(estimates, dims=('stn', 'time'))

                # evaluation
                dtmp1 = stn_value
                dtmp2 = estimates
                metvalue, metname = evaluate_allpoint(dtmp1, dtmp2, 0.1)
                ds_out.coords['met'] = metname
                ds_out[var_poe + '_metric'] = xr.DataArray(metvalue, dims=('stn', 'met'))
                del dtmp1, dtmp2

    # reduce coordinate dims if needed
    if 'x' in ds_out.dims and 'y' in ds_out.dims:
        grid_lat_diff = np.abs(ds_out[grid_lat_name].isel(x=0).values - ds_out[grid_lat_name].isel(x=-1).values)
        grid_lon_diff = np.abs(ds_out[grid_lon_name].isel(y=0).values - ds_out[grid_lon_name].isel(y=-1).values)
        if (np.nanmax(grid_lat_diff) < 1e-10) and (np.nanmax(grid_lon_diff) < 1e-10):
            ds_out.coords['x'] = ds_out[grid_lon_name].isel(y=0).values
            ds_out.coords['y'] = ds_out[grid_lat_name].isel(x=0).values
            ds_out = ds_out.drop_vars([grid_lat_name, grid_lon_name])
            ds_out = ds_out.rename({'y': 'lat', 'x': 'lon'})

    # save output file
    encoding = {}
    for var in ds_out.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}

    ds_out.to_netcdf(outfile, encoding=encoding)

    t2 = time.time()
    print('Time cost (s):', t2-t1)
    print('Regression step completed successfully!\n')

    return config
