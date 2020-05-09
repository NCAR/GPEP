import numpy as np
import netCDF4 as nc
import auxiliary as au
import regression as reg
import datetime as dt
from matplotlib import pyplot as plt
from scipy import io
from scipy.interpolate import interp2d
import os
import sys
from scipy.interpolate import griddata
import h5py
import time
import random


def divide_train_test(data, dividenum, randseed=-1):
    if randseed == -1:
        random.seed(time.time())
    num = len(data)
    subnum = int(num / dividenum)
    data_train = np.zeros([dividenum, num - subnum])
    data_test = np.zeros([dividenum, subnum])
    randindex = random.sample(range(num), num)
    for i in range(dividenum):
        data_test[i, :] = np.sort(data[randindex[i * subnum:(i + 1) * subnum]])
        data_train[i, :] = np.setdiff1d(data, data_test[i])
    return data_train, data_test


def double_cvindex(gmet_stndatafile, dividenum):
    # index for double cross-validation
    datatemp = np.load(gmet_stndatafile)
    prcp_stn0 = datatemp['prcp_stn'][:, 0]
    tmean_stn0 = datatemp['tmean_stn'][:, 0]
    prcp_stnindex = np.argwhere(~np.isnan(prcp_stn0))
    tmean_stnindex = np.argwhere(~np.isnan(tmean_stn0))

    subnum1 = int(len(prcp_stnindex) / dividenum)
    subnum2 = int((len(prcp_stnindex) - subnum1) / dividenum)
    # prcp_testindex1 = np.zeros([dividenum,subnum1])
    # prcp_trainindex1 = np.zeros([dividenum,len(prcp_stnindex) - subnum1])
    prcp_trainindex1, prcp_testindex1 = divide_train_test(prcp_stnindex, dividenum, randseed=123)
    prcp_testindex2 = np.zeros([dividenum, dividenum, subnum2])
    prcp_trainindex2 = np.zeros([dividenum, dividenum, len(prcp_stnindex) - subnum1 - subnum2])
    for i in range(dividenum):
        traini, testi = divide_train_test(prcp_trainindex1[i, :], dividenum, randseed=123)
        prcp_trainindex2[i, :, :] = traini
        prcp_testindex2[i, :, :] = testi

    subnum1 = int(len(tmean_stnindex) / dividenum)
    subnum2 = int((len(tmean_stnindex) - subnum1) / dividenum)
    # tmean_testindex1 = np.zeros([dividenum,subnum1])
    # tmean_trainindex1 = np.zeros([dividenum,len(tmean_stnindex) - subnum1])
    tmean_trainindex1, tmean_testindex1 = divide_train_test(tmean_stnindex, dividenum, randseed=123)
    tmean_testindex2 = np.zeros([dividenum, dividenum, subnum2])
    tmean_trainindex2 = np.zeros([dividenum, dividenum, len(tmean_stnindex) - subnum1 - subnum2])
    for i in range(dividenum):
        traini, testi = divide_train_test(tmean_trainindex1[i, :], dividenum, randseed=123)
        tmean_trainindex2[i, :, :] = traini
        tmean_testindex2[i, :, :] = testi
    return prcp_trainindex1, prcp_testindex1, prcp_trainindex2, prcp_testindex2, \
           tmean_trainindex1, tmean_testindex1, tmean_trainindex2, tmean_testindex2


def calculate_ratio(datatar, dataref, hwsize, upbound=5, lowbound=0.2):
    # datatar, dataref: 2D [nstn, ntime]
    # hwsize: define time window (2*hwsize+1) used to calculate ratio (as ratio for a specific day is too variable)
    # upbound/lowbound: upper and lower limitation of rario
    nearnum = 8  # nearby stations used for extrapolation
    nstn, ntime = np.shape(datatar)

    # 1. calculate ratio for every day
    # lower bound of the ratio
    ratio = np.ones([nstn, ntime])
    for i in range(ntime):
        if i < hwsize:
            windex = np.arange(hwsize * 2 + 1)
        elif i >= ntime - hwsize:
            windex = np.arange(ntime - hwsize * 2 - 1, ntime)
        else:
            windex = np.arange(i - hwsize, i + hwsize + 1)
        dtari = np.nansum(datatar[:, windex], axis=1)
        drefi = np.nansum(dataref[:, windex], axis=1)
        ratio[:, i] = drefi / dtari
    ratio[ratio > upbound] = upbound  # include inf: X / 0
    ratio[ratio < lowbound] = lowbound
    ratio[np.isnan(ratio)] = 1  # 0 / 0
    return ratio


def extrapolation(latin, lonin, datain, latout, lonout, nearnum):
    # datain: one or multiple time steps
    wexp = 3
    if np.ndim(latout) == 1:
        nearstn_loc, nearstn_dist = findnearstn(latin, lonin, latout, lonout, nearnum, 1)
        num = len(latout)
        if np.ndim(datain) == 1:
            datain = datain[:, np.newaxis]
        ntimes = np.shape(datain)[1]
        dataout = np.zeros([num, ntimes])
        for i in range(num):
            dataini = datain[nearstn_loc[i, :], :]
            disti = nearstn_dist[i, :]
            weighti = au.distanceweight(disti, np.max(disti) + 1, wexp)
            weighti = weighti / np.sum(weighti)
            for j in range(ntimes):
                dataout[i, j] = np.sum(dataini[:, j] * weighti)

    elif np.ndim(latout) == 2:
        nearstn_loc, nearstn_dist = findnearstn(latin, lonin, latout, lonout, nearnum, 0)
        if np.ndim(datain) == 2:
            datain = datain[:, :, np.newaxis]
        nrows, ncols, ntimes = np.shape(datain)
        dataout = np.zeros([nrows, ncols, ntimes])
        for r in range(nrows):
            for c in range(ncols):
                dataini = datain[nearstn_loc[r, c, :], :]
                disti = nearstn_dist[r, c, :]
                weighti = au.distanceweight(disti, np.max(disti) + 1, wexp)
                weighti = weighti / np.sum(weighti)
                for j in range(ntimes):
                    dataout[r, c, j] = np.sum(dataini[:, j] * weighti)
    else:
        print('The dimensions of tarlat or tarlon are larger than 2')
        sys.exit()

    return dataout


def findnearstn(stnlat, stnlon, tarlat, tarlon, nearnum, noself):
    # only use lat/lon to find near stations without considering distance in km
    # stnlat/stnlon: 1D
    # tarlat/tarlon: 1D or 2D
    # noself: 1--stnlat and tarlat have overlapped stations, which should be excluded from stnlat

    stnll = np.zeros([len(stnlat), 2])
    stnll[:, 0] = stnlat
    stnll[:, 1] = stnlon

    if len(np.shape(tarlat)) == 1:
        num = len(tarlat)
        nearstn_loc = np.zeros([num, nearnum])
        nearstn_dist = np.zeros([num, nearnum])
        for i in range(num):
            tari = np.array([tarlat[i], tarlon[i]])
            dist = au.distance(tari, stnll)
            if noself == 1:
                dist[dist == 0] = np.inf  # may not be perfect, but work for SCDNA
            indi = np.argsort(dist)
            nearstn_loc[i, :] = indi[0:nearnum]
            nearstn_dist[i, :] = dist[nearstn_loc[i, :]]
    elif len(np.shape(tarlat)) == 2:
        nrows, ncols = np.shape(tarlat)
        nearstn_loc = np.zeros([nrows, ncols, nearnum])
        nearstn_dist = np.zeros([nrows, ncols, nearnum])
        for r in range(nrows):
            for c in range(ncols):
                tari = np.array([tarlat[r, c], tarlon[r, c]])
                dist = au.distance(tari, stnll)
                indi = np.argsort(dist)
                nearstn_loc[r, c, :] = indi[0:nearnum]
                nearstn_dist[r, c, :] = dist[nearstn_loc[r, c, :]]
    else:
        print('The dimensions of tarlat or tarlon are larger than 2')
        sys.exit()

    return nearstn_loc, nearstn_dist


def error_correction(dataori, anomaly, mode='ratio'):
    # default: time is the last dimension
    if mode == 'ratio':
        datacorr = dataori * anomaly
    elif mode == 'diff':
        datacorr = dataori + anomaly
    else:
        sys.exit('Wrong error correction mode')
    return datacorr


def calweight(obsall, preall, mode='RMSE', preprocess=True):
    nstn, ntime = np.shape(obsall)
    met = np.zeros(nstn)
    for i in range(nstn):
        obs = obsall[i, :]
        pre = preall[i, :]
        if preprocess:
            # delete the nan values
            ind_nan = np.isnan(obs) | np.isnan(pre)
            obs = obs[~ind_nan]
            pre = pre[~ind_nan]
        if mode == 'RMSE':
            met[i] = np.sqrt(np.sum(np.square(obs - pre)) / len(obs))  # RMSE
        elif mode == 'CC':
            temp = np.corrcoef(obs, pre)
            met[i] = temp[0][1]  # CC
        else:
            sys.exit('Unknown inputs for calmetric')

    if mode == 'RMSE':
        weight = 1 / (met ** 2)
    elif mode == 'CC':
        met[met < 0] = 0
        weight = met ** 2
    else:
        sys.exit('Unknown inputs for calmetric')

    return weight


def calrmse(dtar, dref):
    if np.ndim(dtar) == 1:
        dtar = dtar[np.newaxis, :]
        dref = dref[np.newaxis, :]
    nstn, ntimes = np.shape(dtar)
    rmse = np.nan * np.zeros(nstn)
    for i in range(nstn):
        rmse[i] = np.sqrt(np.nansum(np.square(dtar[i, :] - dref[i, :])) / ntimes)  # RMSE

    return rmse


def ismember(a, b):
    # tf = np.in1d(a,b) # for newer versions of numpy
    tf = np.array([i in b for i in a])
    u = np.unique(a[tf])
    index = np.array([(np.where(b == i))[0][-1] if t else 0 for i, t in zip(a, tf)])
    return tf, index


def weightmerge(data, weight):
    if np.ndim(data) == 2:
        num, nmodel = np.shape(data)
        dataout = np.zeros(num)
        for i in range(num):
            dataout[i] = np.sum(data[i, :] * weight[i, :]) / np.sum(weight[i, :])
    elif np.ndim(data) == 3:
        nrows, ncols, nmodel = np.shape(data)
        dataout = np.zeros([nrows, ncols])
        for r in range(nrows):
            for c in range(ncols):
                if not np.isnan(data[r, c, 0]):
                    dataout[r, c] = np.sum(data[r, c, :] * weight[r, c, :]) / np.sum(weight[r, c, :])
    return dataout


########################################################################################################################

# basic settings
lontar = np.arange(-180 + 0.05, -50, 0.1)
lattar = np.arange(85 - 0.05, 5, -0.1)
vars = ['prcp', 'tmean', 'trange']
hwsize_ratio = 15  # define time window (2*hwsize+1) used to calculate ratio (as ratio for a specific day is too variable)
nearnum = 8  # the number of nearby stations used to extrapolate points to grids (for correction and merging)
weightmode = 'RMSE'  # the metric used to guide merging (CC or RMSE). Weight = CC**2 or 1/RMSE**2
corrmode = ['ratio', 'diff', 'diff']  # mode for error correction
dividenum = 10  # 10-fold cross-validation

# input files
# station list and data
gmet_stnfile = '/home/gut428/GMET/eCAI_EMDNA/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/home/gut428/stndata_whole.npz'

# downscaled reanalysis data at station points
file_readownstn = ['/ERA5_downto_stn.npz',
                   '/MERRA2_downto_stn.npz',
                   '/JRA55_downto_stn.npz']

# mask file
file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA/DEM/NA_DEM_010deg_trim.mat'

# downscaled reanalysis: gridded data
path_readown = ['', '', '']
prefix = ['ERA5_', 'MERAA2_', 'JRA55_']

# output files
# train and test index file
ttindexfile = '2layer_train_test_index.npz'

# output corrected and merged data
path_reacorr = ['', '', '']
path_merge = ''

########################################################################################################################

# basic processing
# mask
mask = io.loadmat(file_mask)
mask = mask['DEM']
mask[~np.isnan(mask)] = 1  # 1: valid pixels
# attributes
nrows, ncols = np.shape(mask)
lontarm, lattarm = np.meshgrid(lontar, lattar)

########################################################################################################################

# design a two-layer cross-validation: generate station combinations
# index1 extracts 90% stations for merging and 10% stations for validation
# index2 divides the 90% from index1 into 90% and 10% again for error correction
if not os.path.isfile(gmet_stndatafile):
    prcp_trainindex1, prcp_testindex1, prcp_trainindex2, prcp_testindex2, \
    tmean_trainindex1, tmean_testindex1, tmean_trainindex2, tmean_testindex2 = double_cvindex(gmet_stndatafile,
                                                                                              dividenum)
    np.savez_compressed(ttindexfile, prcp_trainindex1=prcp_trainindex1, prcp_testindex1=prcp_testindex1,
                        prcp_trainindex2=prcp_trainindex2, prcp_testindex2=prcp_testindex2,
                        tmean_trainindex1=tmean_trainindex1, tmean_testindex1=tmean_testindex1,
                        tmean_trainindex2=tmean_trainindex2, tmean_testindex2=tmean_testindex2)
    del prcp_trainindex1, prcp_testindex1, prcp_trainindex2, prcp_testindex2, \
        tmean_trainindex1, tmean_testindex1, tmean_trainindex2, tmean_testindex2

taintestindex = np.load(gmet_stndatafile)

########################################################################################################################

# estimate the error of corrected and merged data using cross-validation
# 0. load regression estimates and observations for all stations
# 1. select layer-1 stations (90% train, 10% test)
# 2. select layer-2 stations: 10 combinations (90% train, 10% test) from the 90% stations in layer-1
# 3. for all combinations in layer-2
# 3.1 perform error correction (1) at station points, (2) extrapolate to grids using training stations
# 3.2 perform evaluation using test stations
# 3.3 loop for all combinations and do evaluation at all stations (the 90% training stations from layer-1)
# 3.4 extrapolate evaluation accuracy indicators to the domain
# 4. merge three reanalysis using their indicators in 3.4
# 5. repeat 1-4 for all combinations in layer-1, and get accuracy indicators of the merged dataset for the domain
# 6. use indicators from 3.4 and 5, and select the best one among three reanalysis and merged datasets for each grid
# 7. get the final merged dataset and its accuracy indicators from steo-6

# load regression estimates for all stations
var = 'prcp'  # parameter of modules
reanum = len(file_readownstn)
readata = [''] * reanum
for i in range(reanum):
    di = np.load(file_readownstn[i])
    readata[i] = di[var + '_readown']
    del di

# load observations for all stations
datatemp = np.load(gmet_stndatafile)
stndata = datatemp[var + '_stn']
stnlle = datatemp['stn_lle']
del datatemp

# unify the time length of all data as MERRA2 lacks 1979
stndata = stndata[:, 365:]
for i in range(reanum):
    if np.shape(readata[i])[1] == 14610:  # 1979 to 2018
        readata[i] = readata[i][:, 365:]

nstn, ntimes = np.shape(stndata)
metric_merge_stn = np.nan * np.zeros(nstn)
metric_rea_stn = [''] * reanum
for i in range(reanum):
    metric_rea_stn[i] = np.nan * np.zeros(nstn)

for lay1 in range(dividenum):
    # extract train and test index for layer-1
    if var == 'trange':
        vari = 'tmean'  # trange and tmean have the same index
    else:
        vari = var
    trainindex1 = taintestindex[vari + '_trainindex1'][lay1, :]
    testindex1 = taintestindex[vari + '_testindex1'][lay1, :]
    stndata_trainl1 = stndata[trainindex1, :]
    stndata_testl1 = stndata[testindex1, :]
    stnlle_trainl1 = stnlle[trainindex1, :]
    stnlle_testl1 = stnlle[testindex1, :]

    weight_trainl1 = np.zeros([len(trainindex1), reanum])

    for lay2 in range(dividenum):
        # extract train and test index for layer-2 (subsets of trainindex1)
        trainindex2 = taintestindex[vari + '_trainindex2'][lay1, lay2, :]
        testindex2 = taintestindex[vari + '_testindex2'][lay1, lay2, :]
        stndata_trainl2 = stndata[trainindex2, :]
        stndata_testl2 = stndata[testindex2, :]
        stnlle_trainl2 = stnlle[trainindex2, :]
        stnlle_testl2 = stnlle[testindex2, :]

        for rr in range(reanum):
            readata_trainl2 = readata[rr][trainindex2, :]
            readata_testl2 = readata[rr][testindex2, :]
            # calculate ratio at the train stations
            ratio_ori = calculate_ratio(readata_trainl2, stndata_trainl2, hwsize_ratio, upbound=5, lowbound=0.2)
            # extrapolate the ratio to the test stations (in layer-2)
            ratio_ext = extrapolation(stnlle_trainl2[:, 0], stnlle_trainl2[:, 1], ratio_ori,
                                      stnlle_testl2[:, 0], stnlle_testl2[:, 1], nearnum)
            # correct data at the test stations
            readata_testl2_corr = error_correction(readata_trainl2, ratio_ext, mode=corrmode[rr])
            # estimate weight the test stations
            weight_testl2 = calweight(stndata_testl2, readata_testl2_corr, weightmode)
            # fill the weight to its parent station set (weight_trainl1)
            tf, index = ismember(testindex2, trainindex1)
            weight_trainl1[index, rr] = weight_testl2

    # extrapolate the weight to the test stations (in layer-1)
    weight_testl1 = extrapolation(stnlle_trainl1[:, 0], stnlle_trainl1[:, 1], weight_trainl1,
                                  stnlle_testl1[:, 0], stnlle_testl1[:, 1], nearnum)
    # merge reanalysis products at the test stations
    nstn_testl1 = len(testindex1)
    mergedata_testl1 = np.nan * np.zeros([nstn_testl1, ntimes])
    for i in range(ntimes):
        datain = np.zeros([nstn_testl1, reanum])
        for j in range(reanum):
            datain[:, j] = readata[j][testindex1, i]
        dataout = weightmerge(datain, weight_testl1)
        mergedata_testl1[:, i] = dataout
    # evaluate the merged data at the test stations (just RMSE)
    metric_testl1 = calrmse(mergedata_testl1, stndata_testl1)
    # fill the weight to its parent station set (all original stations)
    metric_merge_stn[testindex1] = metric_testl1

    # repeat error correction using train stations in layer-1 (as in layer-2 only 0.9*0.9=0.81 stations are used)
    for rr in range(reanum):
        readata_trainl1 = readata[rr][trainindex1, :]
        readata_testl1 = readata[rr][testindex1, :]
        ratio_ori = calculate_ratio(readata_trainl1, stndata_trainl1, hwsize_ratio, upbound=5, lowbound=0.2)
        ratio_ext = extrapolation(stnlle_trainl1[:, 0], stnlle_trainl1[:, 1], ratio_ori,
                                  stnlle_testl1[:, 0], stnlle_testl1[:, 1], nearnum)
        readata_testl1_corr = error_correction(readata_trainl1, ratio_ext, mode=corrmode[rr])
        metric_rea_rr = calrmse(readata_testl1_corr, stndata_testl1)
        metric_rea_stn[rr][testindex1] = metric_rea_rr

# extrapolate metric_merge and metric_rea to all grids. the metrics are used to calculate weights
metric_merge_grid = extrapolation(stnlle[:, 0], stnlle[:, 1], metric_merge_stn, lattarm, lontarm, nearnum)

metric_rea_grid = np.zeros([nrows, ncols, reanum])
weight_rea_grid = np.zeros([nrows, ncols, reanum])
for rr in range(reanum):
    metric_rea_grid[:,:,rr] = extrapolation(stnlle[:, 0], stnlle[:, 1], metric_rea_stn[rr], lattarm, lontarm, nearnum)
    weight_rea_grid[:,:,rr] = (1 / (metric_rea_grid[rr] ** 2))
weight_rea_gridsum = np.sum(weight_rea_grid, axis=2)
for rr in range(reanum):
    temp = weight_rea_grid[:,:,rr]/weight_rea_gridsum
    temp[(np.isnan(weight_rea_gridsum)) | (weight_rea_gridsum==0)] = 1/reanum
    weight_rea_grid[:,:,rr] = temp

# which is the best for each grid
merge_choice = np.zeros([nrows, ncols])
weight_all = metric_merge_grid.copy()
for rr in range(reanum):
    weight_all = np.concatenate((weight_all, weight_rea_grid[rr]), axis=2)
merge_choice = np.argmax(weight_all, axis=2) # 0: merge is the best, 1 to N: corresponding corrected reanalysis is the best

########################################################################################################################

# the production of merged product
# load regression estimates for all stations
var = 'prcp'  # parameter of modules
readata = [''] * reanum
for i in range(reanum):
    di = np.load(file_readownstn[i])
    temp = di[var + '_readown']
    if prefix[rr] == 'MERRA2_':
        add = np.nan * np.zeros([nrows, ncols, 365])  # MERRA2 lacks 1979
        temp = np.concatenate((add, temp), axis=2)
    readata[i] = temp
    del di

# load observations for all stations
datatemp = np.load(gmet_stndatafile)
stndata = datatemp[var + '_stn']
stnlle = datatemp['stn_lle']
del datatemp

# error correction of reanalysis based on all stations
dflag = 0
for y in range(1979, 2019):
    readata_gridy = [''] * 3
    fileout_merge = path_merge + '/merge_' + var + '_' + str(y) + '.npz'
    # calculate correction ratio/difference for error correction
    for rr in range(reanum):

        filein = path_readown[rr] + '/' + prefix[rr] + var + '_' + str(y) + '.npz'
        fileout_corr = path_reacorr[rr] + '/' + prefix[rr] + 'corr_' + var + '_' + str(y) + '.npz'

        if prefix[rr] == 'MERRA2_' and y == 1979:
            dayy = 365
            readata_gridy_rr_corr = np.nan * np.zeros([nrows, ncols, dayy])
        else:
            datatemp = np.load(filein)
            readata_gridy_rr = datatemp['data']
            del datatemp
            dayy = np.shape(readata_gridy_rr)[2]
            ratio_ori = calculate_ratio(readata[rr][:, dflag:dflag + dayy], stndata[:, dflag:dflag + dayy],
                                        hwsize_ratio, upbound=5, lowbound=0.2)
            ratio_ext = extrapolation(stnlle[:, 0], stnlle[:, 1], ratio_ori, lattarm, lontarm, nearnum)
            readata_gridy_rr_corr = error_correction(readata_gridy_rr, ratio_ext, mode=corrmode[rr])

        readata_gridy[rr] = readata_gridy_rr_corr
        # save the corrected data
        np.savez_compressed(fileout_corr, datacorr=readata_gridy_rr_corr)

    # merge the reanalysis
    datamerge = np.zeros([nrows,ncols,dayy])
    for rr in range(reanum):
        temp = weight_rea_grid[:,:,rr]
        weightrr = np.repeat(temp[:,:,np.newaxis],dayy,axis=2)
        datamerge = datamerge + weightrr * readata_gridy[rr]
    for r in range(nrows):
        for c in range(ncols):
            if mask[r,c] == 1 and merge_choice[r,c] != 0:
                datamerge[r,c,:] = readata_gridy[merge_choice[r,c]-1][r,c,:]

    np.savez_compressed(fileout_merge,datamerge=datamerge, latitude = lattar, longitude=lontar)
    dflag = dflag + dayy
