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
        data_test[i,:] = np.sort(data[randindex[i * subnum:(i+1) * subnum]])
        data_train[i,:] = np.setdiff1d(data, data_test[i])
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
    return prcp_trainindex1, prcp_testindex1, prcp_trainindex2, prcp_testindex2,\
           tmean_trainindex1, tmean_testindex1, tmean_trainindex2, tmean_testindex2

########################################################################################################################
# station information
# plato
gmet_stnfile = '/home/gut428/GMET/eCAI_EMDNA/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/home/gut428/stndata_whole.npz'
dividenum = 10  # 10-fold cross-validation

########################################################################################################################

# design a two-layer cross-validation: generate station combinations
# index1 extracts 90% stations for merging and 10% stations for validation
# index2 divides the 90% from index1 into 90% and 10% again for error correction
prcp_trainindex1, prcp_testindex1, prcp_trainindex2, prcp_testindex2, \
tmean_trainindex1, tmean_testindex1, tmean_trainindex2, tmean_testindex2 = double_cvindex(gmet_stndatafile, dividenum)

########################################################################################################################
# load regression estimates for all stations
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



