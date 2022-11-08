# correct the bias in reanalysis products and merge different reanalysis products
# No correction. Just merge based on downscaled reanalysis data


# import numpy as np
# import auxiliary as au
# from matplotlib import pyplot as plt
from scipy import io
import os
import netCDF4 as nc
from bma_merge import bma
from auxiliary_merge import *


def empirical_cdf(data, probtar):
    # data: vector of data
    data2 = data[~np.isnan(data)]
    if len(data2) > 0:
        ds = np.sort(data2)
        probreal = np.arange(len(data2)) / (len(data2) + 1)
        ecdf_out = np.interp(probtar, probreal, ds)
    else:
        ecdf_out = np.nan * np.zeros(len(probtar))
    return ecdf_out


def cdf_correction(cdf_ref, value_ref, cdf_raw, value_raw, value_tar):
    prob_tar = np.interp(value_tar, value_raw, cdf_raw)
    value_out = np.interp(prob_tar, cdf_ref, value_ref)
    return value_out


def calculate_anomaly(datatar, dataref, hwsize, amode, upbound=5, lowbound=0.2):
    # datatar, dataref: 2D [nstn, ntime]
    # amode: anomaly mode ('ratio' or 'diff')
    # hwsize: define time window (2*hwsize+1) used to calculate ratio (as ratio for a specific day is too variable)
    # upbound/lowbound: upper and lower limitation of ratio/difference
    if np.ndim(datatar) == 1:  # only one time step
        datatar = datatar[:, np.newaxis]
        dataref = dataref[:, np.newaxis]

    nstn, ntime = np.shape(datatar)
    if ntime < hwsize * 2 + 1:
        print('The window size is larger than time steps when calculating ratio between tar and ref datasets')
        print('Please set a smaller hwsize')
        sys.exit()

    anom = np.ones([nstn, ntime])

    for i in range(ntime):
        if i < hwsize:
            windex = np.arange(hwsize * 2 + 1)
        elif i >= ntime - hwsize:
            windex = np.arange(ntime - hwsize * 2 - 1, ntime)
        else:
            windex = np.arange(i - hwsize, i + hwsize + 1)
        dtari = np.nanmean(datatar[:, windex], axis=1)
        drefi = np.nanmean(dataref[:, windex], axis=1)

        if amode == 'ratio':
            temp = drefi / dtari
            temp[(dtari == 0) & (drefi == 0)] = 1
            anom[:, i] = temp
        elif amode == 'diff':
            anom[:, i] = drefi - dtari
        else:
            sys.exit('Unknow amode. Please use either ratio or diff')

    anom[anom > upbound] = upbound
    anom[anom < lowbound] = lowbound
    return anom


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
        nearstn_loc = -1 * np.ones([num, nearnum], dtype=int)
        nearstn_dist = -1 * np.ones([num, nearnum], dtype=float)
        for i in range(num):
            if np.mod(i, 100) == 0:
                print('station', i, nstn)
            if np.isnan(tarlat[i]) or np.isnan(tarlon[i]):
                continue
            tari = np.array([tarlat[i], tarlon[i]])
            dist = au.distance(tari, stnll)
            dist[np.isnan(dist)] = 1000000000
            if noself == 1:
                dist[dist == 0] = np.inf  # may not be perfect, but work for SCDNA
            indi = np.argsort(dist)
            nearstn_loc[i, :] = indi[0:nearnum]
            nearstn_dist[i, :] = dist[nearstn_loc[i, :]]
    elif len(np.shape(tarlat)) == 2:
        nrows, ncols = np.shape(tarlat)
        nearstn_loc = -1 * np.ones([nrows, ncols, nearnum], dtype=int)
        nearstn_dist = -1 * np.ones([nrows, ncols, nearnum], dtype=float)
        for r in range(nrows):
            print('rows', r, nrows)
            for c in range(ncols):
                if np.isnan(tarlat[r, c]) or np.isnan(tarlon[r, c]):
                    continue
                tari = np.array([tarlat[r, c], tarlon[r, c]])
                dist = au.distance(tari, stnll)
                dist[np.isnan(dist)] = 1000000000
                indi = np.argsort(dist)
                nearstn_loc[r, c, :] = indi[0:nearnum]
                nearstn_dist[r, c, :] = dist[nearstn_loc[r, c, :]]
    else:
        print('The dimensions of tarlat or tarlon are larger than 2')
        sys.exit()

    return nearstn_loc, nearstn_dist


def error_correction_stn(corrmode, stndata_i2_near, nearstn_weighti2, readata_stn_i2, readata_i2_near, ecdf_prob):
    # corrmode: QM, Mul_Climo, Mul_Daily, Add_Climo, Add_Climo
    nearstn_numi2, ntimes = np.shape(stndata_i2_near)
    corrdata_out = np.zeros(ntimes)
    if corrmode == 'QM':
        cdf_rea = empirical_cdf(readata_stn_i2, ecdf_prob)
        for j in range(nearstn_numi2):
            cdf_ref = empirical_cdf(stndata_i2_near[j, :], ecdf_prob)
            qmdata_rj = cdf_correction(ecdf_prob, cdf_ref, ecdf_prob, cdf_rea, readata_stn_i2)
            corrdata_out = corrdata_out + qmdata_rj * nearstn_weighti2[j]
        corrdata_out = corrdata_out / np.sum(nearstn_weighti2)
    elif corrmode[0:3] == 'Mul' or corrmode[0:3] == 'Add':
        # multplicative correction or additive correction
        if corrmode[4:] == 'Daily':
            dtar = readata_i2_near
            dref = stndata_i2_near
        elif corrmode[4:] == 'Climo':
            dtar = np.nanmean(readata_i2_near, axis=1)
            dtar = dtar[:, np.newaxis]
            dref = np.nanmean(stndata_i2_near, axis=1)
            dref = dref[:, np.newaxis]
        else:
            sys.exit('Unknown corrmode')
        if corrmode[0:3] == 'Mul':
            corrfactor_i2_near = calculate_anomaly(dtar, dref, 0, 'ratio', 10, 0)  # 10 is default max limit
        else:
            corrfactor_i2_near = calculate_anomaly(dtar, dref, 0, 'diff', 9999, -9999)
        weight_use = np.tile(nearstn_weighti2, (np.shape(corrfactor_i2_near)[1], 1)).T
        weight_use[np.isnan(corrfactor_i2_near)] = np.nan
        corrfactor_i2 = np.nansum(corrfactor_i2_near * weight_use, axis=0) / np.nansum(weight_use)
        if corrmode[0:3] == 'Mul':
            corrdata_out = readata_stn_i2 * corrfactor_i2
        else:
            corrdata_out = readata_stn_i2 + corrfactor_i2
    else:
        sys.exit('Unknown corrmode')

    return corrdata_out


def error_correction(dataori, anomaly, mode='ratio'):
    # default: time is the last dimension
    if mode == 'ratio':
        datacorr = dataori * anomaly
    elif mode == 'diff':
        datacorr = dataori + anomaly
    else:
        sys.exit('Wrong error correction mode')
    return datacorr


def calweight(obs, rea, mode, preprocess=True):
    ntimes, reanum = np.shape(rea)
    if preprocess:
        # delete the nan values
        ind_nan = np.isnan(obs + np.sum(rea, axis=1))
        obs = obs[~ind_nan]
        rea = rea[~ind_nan, :]

    if len(obs) > 2:
        if mode == 'BMA':
            weight, sigma, sigma_s = bma(rea, obs)
        else:
            met = np.zeros(reanum)
            if mode == 'RMSE':
                for i in range(reanum):
                    met[i] = np.sqrt(np.sum(np.square(obs - rea[:, i])) / len(obs))  # RMSE
                weight = 1 / (met ** 2)
            elif mode == 'CC':
                for i in range(reanum):
                    met[i] = np.corrcoef(obs, rea[:, i])[0][1]
                weight = (met ** 2)
    else:
        weight = np.ones(reanum) / reanum
    if np.any(np.isnan(weight)):
        weight = np.ones(reanum) / reanum
    return weight


def weightmerge(data, weight):
    if np.ndim(data) == 2:
        weight2 = weight.copy()
        weight2[np.isnan(data)] = np.nan
        dataout = np.nansum(data * weight2, axis=1) / np.nansum(weight2, axis=1)
    elif np.ndim(data) == 3:
        weight2 = weight.copy()
        weight2[np.isnan(data)] = np.nan
        dataout = np.nansum(data * weight2, axis=2) / np.nansum(weight2, axis=2)
        dataout[np.isnan(data[:, :, 0])] = np.nan
    return dataout



def merge_stn(stndata, readata_stn, nearstn_loc, nearstn_dist, var, weightmode):
    reanum, nstn, ntimes = np.shape(readata_stn)

    # initialization
    reacorr_stn = readata_stn  # corrected reanalysis data
    reamerge_weight_stn = np.nan * np.zeros([nstn, reanum], dtype=np.float32)  # weight used to obtain reamerge_stn
    reamerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)  # merged reanalysis at station points

    # get merge weight for all stations
    print('calculate merging weights')
    for i1 in range(nstn):
        if np.mod(i1, 1000) == 0:
            print(i1)
        stndata_i1 = stndata[i1, :]
        corrdata_i1 = reacorr_stn[:, i1, :].T
        # get the final merging weight
        if weightmode == 'BMA' and var == 'prcp':
            # exclude zero precipitation and carry out box-cox transformation
            datatemp = np.zeros([ntimes, reanum + 1])
            datatemp[:, 0] = stndata_i1
            datatemp[:, 1:] = corrdata_i1
            ind0 = np.sum(datatemp >= 0.01, axis=1) == (reanum + 1)  # positive hit events
            dobs = box_cox_transform(stndata_i1[ind0])
            drea = box_cox_transform(corrdata_i1[ind0, :])
        else:
            dobs = stndata_i1
            drea = corrdata_i1
        reamerge_weight_stn[i1, :] = calweight(dobs, drea, weightmode)

    print('merge data using leave-one-out')
    for i1 in range(nstn):  # layer-1
        if np.mod(i1, 1000) == 0:
            print(i1)
        if np.isnan(stndata[i1, 0]):
            continue
        nearstn_loci1 = nearstn_loc[i1, :]
        nearstn_disti1 = nearstn_dist[i1, :]
        induse = nearstn_loci1 > -1
        nearstn_loci1 = nearstn_loci1[induse]
        nearstn_disti1 = nearstn_disti1[induse]
        nearstn_numi1 = len(nearstn_loci1)
        if nearstn_numi1 == 0:
            sys.exit('No near station for the target station (layer-1)')

        maxd = np.max([np.max(nearstn_disti1) + 1, 100])
        nearstn_weighti1 = au.distanceweight(nearstn_disti1, maxd, 3)
        nearstn_weighti1 = nearstn_weighti1 / np.sum(nearstn_weighti1)

        reamerge_weight_i2 = reamerge_weight_stn[nearstn_loci1, :]
        # get merging weight at i1 and merge reanalysis
        # note: this weight is just for independent merging so we can estimate the error of merged reanalysis
        # the real weight will be estimated using just one-layer cross-validation
        weight_use = np.tile(nearstn_weighti1, (reanum, 1)).T
        weight_i1 = np.sum(weight_use * reamerge_weight_i2, axis=0)
        weight_i1 = weight_i1 / np.sum(weight_i1)

        weight_use = np.tile(weight_i1, (ntimes, 1))
        weight_use[np.isnan(corrdata_i1)] = np.nan
        reamerge_stni1 = np.nansum(weight_use * corrdata_i1, axis=1) / np.nansum(weight_use, axis=1)
        reamerge_stn[i1, :] = reamerge_stni1

    # note: reamerge_weight_stn is the final merging weight, and reacorr_stn is the final corrected data
    # but reamerge_stn is just independent merging estimates which is calculated from 2-layer cross validation
    return reamerge_stn, reamerge_weight_stn, reacorr_stn


def correction_rea(stndata, ecdf_prob, readata_stn, nearstn_loc, nearstn_dist, corrmode):
    # compare the performance of daily-scale multiplicative correction and QM
    reanum, nstn, ntimes = np.shape(readata_stn)
    nprob = len(ecdf_prob)

    # initialization
    reacorr = np.nan * np.zeros([reanum, nstn, ntimes], dtype=np.float32)  # corrected reanalysis data

    for i1 in range(nstn):  # layer-1
        # if np.mod(i1, 1000) == 0:
        #     print(i1)
        if np.isnan(stndata[i1, 0]):
            continue

        nearstn_loci1 = nearstn_loc[i1, :]
        nearstn_disti1 = nearstn_dist[i1, :]
        induse = nearstn_loci1 > -1
        nearstn_loci1 = nearstn_loci1[induse]
        nearstn_disti1 = nearstn_disti1[induse]
        nearstn_numi1 = len(nearstn_loci1)
        if nearstn_numi1 == 0:
            sys.exit('No near station for the target station (layer-1)')

        stndata_i1_near = stndata[nearstn_loci1, :]
        readata_stn_i1 = readata_stn[:, i1, :]
        readata_i1_near = readata_stn[:, nearstn_loci1, :]
        maxd = np.max([np.max(nearstn_disti1) + 1, 100])
        nearstn_weighti1 = au.distanceweight(nearstn_disti1, maxd, 3)
        nearstn_weighti1 = nearstn_weighti1 / np.sum(nearstn_weighti1)

        # get corrected data at i1
        corrdata_i1 = np.zeros([ntimes, reanum])
        for r in range(reanum):
            corrdata_i1[:, r] = error_correction_stn(corrmode, stndata_i1_near, nearstn_weighti1,
                                                     readata_stn_i1[r, :], readata_i1_near[r, :, :], ecdf_prob)
        reacorr[:, i1, :] = corrdata_i1.T

    return reacorr


########################################################################################################################
# time periods and methods
# var = 'prcp'  # ['prcp', 'tmean', 'trange']: this should be input from sbtach script
# weightmode for merging: (CC, RMSE, BMA): Weight = CC**2, or Weight = 1/RMSE**2, or Weight = BMA
# corrmode: QM, Mul_Climo, Mul_Daily, Add_Climo, Add_Climo
# year range for merging. note weight is calculated using all data not limited by year

# read from inputs
var = sys.argv[1]
weightmode = sys.argv[2]
corrmode = sys.argv[3]
y1 = int(sys.argv[4])
y2 = int(sys.argv[5])
year = [y1, y2] # this is only useful for gridded correction/merge

# var = sys.argv[1]
# weightmode = sys.argv[2]
# corrmode = sys.argv[3]
# month = int(sys.argv[4])

# embeded
# var = 'prcp'
# weightmode = 'BMA'
# corrmode = 'QM'
# y1 = 2017
# y2 = 2018
# year = [y1, y2]

print('var is ', var)
print('weightmode is ', weightmode)
print('years are ', y1, y2)

########################################################################################################################

# basic settings
nearnum = 10  # the number of nearby stations used to extrapolate points to grids (for correction and merging)
prefix = ['ERA5_', 'MERRA2_', 'JRA55_']

### Plato settings
# input files/paths
FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/gridinfo_whole.nc'
gmet_stnfile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_aftercheck.npz'
file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA_new/DEM/NA_DEM_010deg_trim.mat'
path_readowngrid = ['/datastore/GLOBALWATER/CommonData/EMDNA_new/ERA5_day_ds',
                    '/datastore/GLOBALWATER/CommonData/EMDNA_new/MERRA2_day_ds',
                    '/datastore/GLOBALWATER/CommonData/EMDNA_new/JRA55_day_ds']
file_readownstn = ['/datastore/GLOBALWATER/CommonData/EMDNA_new/ERA5_day_ds/ERA5_downto_stn_GWR.npz',
                   '/datastore/GLOBALWATER/CommonData/EMDNA_new/MERRA2_day_ds/MERRA2_downto_stn_GWR.npz',
                   '/datastore/GLOBALWATER/CommonData/EMDNA_new/JRA55_day_ds/JRA55_downto_stn_GWR.npz']
file_nearstn = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck/nearstn_catalog.npz'

# output files/paths (can also be used as inputs once generated)
path_merge = '/home/gut428/ReanalysisCorrMerge/GWRBMA_merge'
### Plato settings

file_corrmerge_stn = path_merge + '/mergecorr_stn_' + var + '_GWR_' + weightmode + '.npz'  # file of indepedent corrected/merging data and merging weights

########################################################################################################################

# basic processing
print('start basic processing')
reanum = len(file_readownstn)

# mask
mask = io.loadmat(file_mask)
mask = mask['DEM']
mask[~np.isnan(mask)] = 1  # 1: valid pixels
nrows, ncols = np.shape(mask)

# meshed lat/lon of the target region
ncfid = nc.Dataset(FileGridInfo)
lattarm = ncfid.variables['latitude'][:].data
lattarm = np.flipud(lattarm)
lontarm = ncfid.variables['longitude'][:].data
ncfid.close()
lontarm[np.isnan(mask)] = np.nan
lattarm[np.isnan(mask)] = np.nan
lontar = lontarm[0, :]
lattar = lattarm[:, 0]

# load observations for all stations
datatemp = np.load(gmet_stndatafile)
stndata = datatemp[var + '_stn']
stninfo = datatemp['stninfo']
stnlle = stninfo[:,1:4]
date_ymd = datatemp['date_ymd']
nstn, ntimes = np.shape(stndata)
del datatemp
date_yyyy = (date_ymd/10000).astype(int)
date_mm = (np.mod(date_ymd, 10000)/100).astype(int)

# probability bins for QM
binprob = 500
ecdf_prob = np.arange(0, 1 + 1 / binprob, 1 / binprob)

########################################################################################################################

# load near station information
if not os.path.isfile(file_nearstn):
    print(file_nearstn, 'does not exist')
    sys.exit()
datatemp = np.load(file_nearstn)
if var == 'prcp':
    nearstn_loc = datatemp['near_stn_prcpLoc'][:, 0:nearnum]
    nearstn_dist = datatemp['near_stn_prcpDist'][:, 0:nearnum]
    neargrid_loc = datatemp['near_grid_prcpLoc'][:, :, 0:nearnum]
    neargrid_dist = datatemp['near_grid_prcpDist'][:, :, 0:nearnum]
else:
    nearstn_loc = datatemp['near_stn_tempLoc'][:, 0:nearnum]
    nearstn_dist = datatemp['near_stn_tempDist'][:, 0:nearnum]
    neargrid_loc = datatemp['near_grid_tempLoc'][:, :, 0:nearnum]
    neargrid_dist = datatemp['near_grid_tempDist'][:, :, 0:nearnum]

# this is because of the inverse latitude rule
neargrid_loc = np.flipud(neargrid_loc)
neargrid_dist = np.flipud(neargrid_dist)

########################################################################################################################

# load downscaled reanalysis at station points
print('load downscaled reanalysis data at station points')
readata_stn = np.nan * np.zeros([reanum, nstn, ntimes], dtype=np.float32)
for rr in range(reanum):
    dr = np.load(file_readownstn[rr])
    temp = dr[var + '_readown']
    readata_stn[rr, :, :] = temp
    del dr, temp
if var == 'prcp':
    readata_stn[readata_stn < 0] = 0

########################################################################################################################

# get merged reanalysis data at all station points using two-layer cross-validation
if os.path.isfile(file_corrmerge_stn):
    print('load independent merged/corrected data at station points')
    datatemp = np.load(file_corrmerge_stn)
    reamerge_stn = datatemp['reamerge_stn']
    reamerge_weight_stn = datatemp['reamerge_weight_stn']
    reacorr_stn = datatemp['reacorr_stn']
    del datatemp
else:
    print('calculate independent merged/corrected data at station points')
    reamerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)
    reamerge_weight_stn = np.nan * np.zeros([12, nstn, reanum], dtype=np.float32)
    reacorr_stn = np.nan * np.zeros([reanum, nstn, ntimes], dtype=np.float32)
    # for each month
    for m in range(12):
        print('month', m + 1)
        indm = date_mm == (m + 1)
        reamerge_stnm, reamerge_weight_stnm, reacorr_stnm = \
            merge_stn(stndata[:, indm], readata_stn[:, :, indm], nearstn_loc, nearstn_dist, var, weightmode)
        reamerge_stn[:, indm] = reamerge_stnm
        reacorr_stn[:, :, indm] = reacorr_stnm
        reamerge_weight_stn[m, :, :] = reamerge_weight_stnm

    # the variables are independent with their concurrent stations. thus, station data can be used to evaluate them
    np.savez_compressed(file_corrmerge_stn, reamerge_stn=reamerge_stn, reamerge_weight_stn=reamerge_weight_stn,
                        reacorr_stn=reacorr_stn, date_ymd=date_ymd, prefix=prefix, stninfo=stninfo)

########################################################################################################################

# interpolate merging weights to grids
filebma_merge_weight = path_merge + '/mergeweight_' + var + '_' + weightmode + '.npz'
if os.path.isfile(filebma_merge_weight):
    print('Load merging weight')
    datatemp = np.load(filebma_merge_weight)
    reamerge_weight_grid = datatemp['reamerge_weight_grid']
    del datatemp
else:
    print('Interpolate merging weight')
    reamerge_weight_grid = np.nan * np.zeros([12, reanum, nrows, ncols], dtype=np.float32)
    for m in range(12):
        for rr in range(reanum):
            reamerge_weight_grid[m, rr, :, :] = extrapolation(reamerge_weight_stn[m, :, rr], neargrid_loc,
                                                              neargrid_dist)
    np.savez_compressed(filebma_merge_weight, reamerge_weight_grid=reamerge_weight_grid, reaname=prefix,
                        latitude=lattar, longitude=lontar)

#######################################################################################################################


# start final correction and merging

# process for each month
for m in range(12):
    print('Correction and Merge: month', m + 1)
    # correction and merging for each year
    for y in range(year[0], year[1] + 1):
        print('Year', y)
        filebma_merge = path_merge + '/bmamerge_' + var + '_' + str(y * 100 + m + 1) + '.npz'
        if os.path.isfile(filebma_merge):
            print('file exists ... continue')
            continue
        # date processing
        mmy = date_mm.copy()
        mmy = mmy[date_yyyy == y]
        indmmy = mmy == (m + 1)
        indmmy2 = (date_yyyy == y) & (date_mm == m + 1)
        mmdays = np.sum(indmmy)

        # read raw gridded reanalysis data
        readata_raw = np.nan * np.zeros([reanum, nrows, ncols, mmdays], dtype=np.float32)
        for rr in range(reanum):
            if not (prefix[rr] == 'MERRA2_' and y == 1979):
                filer = path_readowngrid[rr] + '/' + prefix[rr] + 'ds_' + var + '_' + str(y * 100 + m + 1) + '.npz'
                d = np.load(filer)
                readata_raw[rr, :, :, :] = d['data']
                del d

        ################################################################################################################
        print('Reanalysis merging')
        # start BMA-based merging
        if not os.path.isfile(filebma_merge):
            # initialization
            bma_data = np.nan * np.zeros([nrows, ncols, mmdays], dtype=np.float32)
            # bma_error = np.nan * np.zeros([nrows, ncols, mmdays], dtype=np.float32)

            # (1) estimate the error of corrected data by interpolating stations
            bma_error = extrapolation(reamerge_stn[:, indmmy2] - stndata[:, indmmy2], neargrid_loc, neargrid_dist)

            # (2) estimate the value of merged data
            reamerge_weight_gridm = reamerge_weight_grid[m, :, :, :].copy()
            for i in range(mmdays):
                datai = readata_raw[:, :, :, i]
                weighti = reamerge_weight_gridm.copy()
                weighti[np.isnan(datai)] = np.nan
                bma_data[:, :, i] = np.nansum(weighti * datai, axis=0) / np.nansum(weighti, axis=0)
            np.savez_compressed(filebma_merge, bma_data=bma_data, bma_error=bma_error,
                                reaname=prefix, latitude=lattar, longitude=lontar)

            del bma_error, bma_data