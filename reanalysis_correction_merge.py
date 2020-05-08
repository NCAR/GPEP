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
    wexp = 3
    if np.ndim(latout) == 1:
        nearstn_loc, nearstn_dist = findnearstn(latin, lonin, latout, lonout, nearnum, 1)
        num = len(latout)
        if np.ndim(datain) == 1:
            datain = datain[:,np.newaxis]
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
            datain = datain[:,:,np.newaxis]
        nrows, ncols, ntimes = np.shape(datain)
        dataout = np.zeros([nrows, ncols, ntimes])
        for r in range(nrows):
            for c in range(ncols):
                dataini = datain[nearstn_loc[r,c, :], :]
                disti = nearstn_dist[r,c,:]
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

def error_correction(dataori, anomaly, mode = 'ratio'):
    # default: time is the last dimension
    if mode == 'ratio':
        datacorr = dataori * anomaly
    elif mode == 'diff':
        datacorr = dataori + anomaly
    else:
        sys.exit('Wrong error correction mode')
    return datacorr

def calweight(obsall, preall, mode = 'RMSE', preprocess=True):
    nstn, ntime = np.shape(obsall)
    met = np.zeros(nstn)
    for i in range(nstn):
        obs = obsall[i,:]
        pre = preall[i,:]
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
        weight = 1 / (met**2)
    elif mode == 'CC':
        met[met<0] = 0
        weight = met ** 2
    else:
        sys.exit('Unknown inputs for calmetric')

    return weight

def ismember(a, b):
    # tf = np.in1d(a,b) # for newer versions of numpy
    tf = np.array([i in b for i in a])
    u = np.unique(a[tf])
    index = np.array([(np.where(b == i))[0][-1] if t else 0 for i,t in zip(a,tf)])
    return tf, index


def weightmerge(data,weight):
    pass


########################################################################################################################

# basic settings
lontar = np.arange(-180 + 0.05, -50, 0.1)
lattar = np.arange(85 - 0.05, 5, -0.1)
vars = ['prcp', 'tmean', 'trange']
hwsize_ratio = 15 # define time window (2*hwsize+1) used to calculate ratio (as ratio for a specific day is too variable)
nearnum = 8 # the number of nearby stations used to extrapolate points to grids (for correction and merging)
weightmode = 'RMSE' # the metric used to guide merging (CC or RMSE). Weight = CC**2 or 1/RMSE**2

########################################################################################################################

# station information
gmet_stnfile = '/home/gut428/GMET/eCAI_EMDNA/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/home/gut428/stndata_whole.npz'

# train and test index file
ttindexfile = '2layer_train_test_index.npz'
dividenum = 10  # 10-fold cross-validation

# downscaled reanalysis data at station points
file_readownstn = ['/ERA5_downto_stn.npz',
                   '/MERRA2_downto_stn.npz',
                   '/JRA55_downto_stn.npz']
reanum = len(file_readownstn)

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

# 0. load regression estimates for all stations
var = 'prcp'  # parameter of modules
readata = [''] * reanum
for i in range(reanum):
    di = np.load(file_readownstn[i])
    readata[i] = di[var + '_readown']
    del di

# 0. load observations for all stations
datatemp = np.load(gmet_stndatafile)
stndata = datatemp[var + '_stn']
stnlle = datatemp['stn_lle']
del datatemp

for lay1 in range(dividenum):
    # extract train and test index for layer-1
    if var == 'trange':
        vari = 'tmean'  # trange and tmean have the same index
    else:
        vari = var
    trainindex1 = taintestindex[vari + '_trainindex1'][lay1, :]
    testindex1 = taintestindex[vari + '_testindex1'][lay1, :]
    stnlle_trainl1 = stnlle[trainindex1, :]
    stnlle_testl1 = stnlle[testindex1, :]

    weight_trainl1 = np.zeros([len(trainindex1),reanum])

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
            ratio_ext = extrapolation(stnlle_trainl2[:,0], stnlle_trainl2[:,1], ratio_ori,
                                       stnlle_testl2[:,0], stnlle_testl2[:,1], nearnum)
            # correct data at the test stations
            readata_testl2_corr = error_correction(readata_trainl2, ratio_ext, mode = 'ratio')
            # estimate weight the test stations
            weight_testl2  = calweight(stndata_testl2, readata_testl2_corr, weightmode)
            # fill the weight to its parent station set (weight_trainl1)
            tf, index = ismember(testindex2, trainindex1)
            weight_trainl1[index,rr] = weight_testl2

    # extrapolate the weight to the test stations (in layer-1)
    weight_testl1 = extrapolation(stnlle_trainl1[:, 0], stnlle_trainl1[:, 1], weight_trainl1,
                              stnlle_testl1[:, 0], stnlle_testl1[:, 1], nearnum)
    # merge reanalysis products at the test stations


