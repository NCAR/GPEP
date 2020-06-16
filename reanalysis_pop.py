import numpy as np
import regression as reg
from scipy import io
from auxiliary_merge import m_DateList
import os, sys
from bma_merge import bma
from auxiliary_merge import extrapolation

# other choices for logistic regression
# LogisticRegression with solver='lbfgs' is two times faster than the least-square iterations
# SGDClassifier with default setting is sevent times faster, but not as accurate
# more testing would be needed
from sklearn.linear_model import LogisticRegression, SGDClassifier


# read from inputs
# time1 = int(sys.argv[1])
# time2 = int(sys.argv[2])
# print(time1,time2)

yearin = int(sys.argv[1])
monthin = int(sys.argv[2])
print(yearin,monthin)


prefix = ['ERA5_', 'MERRA2_', 'JRA55_']

# ### Local Mac settings
# # input files/paths
# gmet_stnfile = '/Users/localuser/Research/EMDNA/basicinfo/stnlist_whole.txt'
# gmet_stndatafile = '/Users/localuser/Research/EMDNA/stndata_whole.npz'
# file_mask = './DEM/NA_DEM_010deg_trim.mat'
# near_file_GMET = '/Users/localuser/Research/EMDNA/regression/weight_nearstn.npz' # near station of stations/grids
# path_readowngrid = ['/Users/localuser/Research/EMDNA/downscale/ERA5',  # downscaled gridded data
#                     '/Users/localuser/Research/EMDNA/downscale/MERRA2',
#                     '/Users/localuser/Research/EMDNA/downscale/JRA55']
# file_readownstn = ['/Users/localuser/Research/EMDNA/downscale/ERA5_downto_stn_nearest.npz', # downscaled to stn points
#                    '/Users/localuser/Research/EMDNA/downscale/MERRA2_downto_stn_nearest.npz',
#                    '/Users/localuser/Research/EMDNA/downscale/JRA55_downto_stn_nearest.npz']
#
# # output files/paths (can also be used as inputs once generated)
# path_ecdf = '/Users/localuser/Research/EMDNA/merge/ECDF'
# path_pop = '/Users/localuser/Research/EMDNA/pop'
# ### Local Mac settings


### Plato settings
gmet_stnfile = '/datastore/GLOBALWATER/CommonData/EMDNA/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA/stndata_whole.npz'
file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA/DEM/NA_DEM_010deg_trim.mat'
near_file_GMET = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout/weight.npz'
path_readowngrid = ['/datastore/GLOBALWATER/CommonData/EMDNA/ERA5_day_ds',  # downscaled gridded data
                   '/datastore/GLOBALWATER/CommonData/EMDNA/MERRA2_day_ds',
                   '/datastore/GLOBALWATER/CommonData/EMDNA/JRA55_day_ds']
file_readownstn = ['/datastore/GLOBALWATER/CommonData/EMDNA/ERA5_day_ds/ERA5_downto_stn_GWR.npz', # downscaled to stn points
                   '/datastore/GLOBALWATER/CommonData/EMDNA/MERRA2_day_ds/MERRA2_downto_stn_GWR.npz',
                   '/datastore/GLOBALWATER/CommonData/EMDNA/JRA55_day_ds/JRA55_downto_stn_GWR.npz']
path_ecdf = '/datastore/GLOBALWATER/CommonData/EMDNA/ReanalysisCorrMerge/ECDF'
path_pop = '/home/gut428/ReanalysisCorrMerge/pop'
### Plato settings

file_reapop_stn = path_pop + '/reanalysis_pop_stn.npz'
file_popmerge_stn = path_pop + '/bmamerge_pop_stn.npz'

########################################################################################################################

# basic processing
print('start basic processing')

lontar = np.arange(-180 + 0.05, -50, 0.1)
lattar = np.arange(85 - 0.05, 5, -0.1)
# mask
mask = io.loadmat(file_mask)
mask = mask['DEM']
mask[~np.isnan(mask)] = 1  # 1: valid pixels

# meshed lat/lon of the target region
reanum = len(file_readownstn)
nrows, ncols = np.shape(mask)
lontarm, lattarm = np.meshgrid(lontar, lattar)
lontarm[np.isnan(mask)] = np.nan
lattarm[np.isnan(mask)] = np.nan

# date list
date_list, date_number = m_DateList(1979, 2018, 'ByYear')

# load observations for all stations
datatemp = np.load(gmet_stndatafile)
stndata = datatemp['prcp_stn']
stnlle = datatemp['stn_lle']
nstn, ntimes = np.shape(stndata)
del datatemp

# load near station information
datatemp = np.load(near_file_GMET)
near_loc_stn = datatemp['near_stn_prcpLoc']
near_weight_stn = datatemp['near_stn_prcpWeight']
near_dist_stn = datatemp['near_stn_prcpDist']
near_loc_grid = datatemp['near_grid_prcpLoc']
near_weight_grid = datatemp['near_grid_prcpWeight']
near_dist_grid = datatemp['near_grid_prcpDist']
near_loc_grid = np.flipud(near_loc_grid)
near_weight_grid = np.flipud(near_weight_grid)
near_dist_grid = np.flipud(near_dist_grid)

# probability bins for QM
binprob = 500
ecdf_prob = np.arange(0, 1 + 1 / binprob, 1 / binprob)

########################################################################################################################

# load downscaled reanalysis at station points
print('load downscaled reanalysis data at station points')
readata_stn = np.nan * np.zeros([reanum, nstn, ntimes], dtype=np.float32)
for rr in range(reanum):
    dr = np.load(file_readownstn[rr])
    temp = dr['prcp_readown']
    readata_stn[rr, :, :] = temp
    del dr, temp
readata_stn[readata_stn < 0] = 0

########################################################################################################################

# method-1: estimate pop using a univariate regression between station occurrence (0-1) and reanalysis precipitation
# file_popt = path_pop + '/reapop_stn_' + str(time1) + '-' + str(time2) + '.npz'
if os.path.isfile(file_reapop_stn):
    datatemp = np.load(file_reapop_stn)
    reapop_stn = datatemp['reapop_stn']
    del datatemp
else:
    reapop_stn = np.nan * np.zeros([reanum, nstn, ntimes], dtype=np.float32)
    for rr in range(reanum):
        print('reanalysis',rr)
        for gg in range(nstn):
            if np.mod(gg,100)==0:
                print(gg,nstn)
            if np.isnan(stndata[gg, 0]):
                continue
            nearloc = near_loc_stn[gg, :]
            neardist = near_dist_stn[gg, :]
            nearweight = near_weight_stn[gg, :]
            neardist = neardist[nearloc > -1]
            nearweight = nearweight[nearloc > -1]
            nearweight = nearweight/np.sum(nearweight)
            nearloc = nearloc[nearloc > -1]

            nstn_prcp = len(nearloc)
            w_pcp_red = np.zeros([nstn_prcp, nstn_prcp])
            for i in range(nstn_prcp):
                w_pcp_red[i, i] = nearweight[i]  # eye matrix: stn weight in one-one lien

            x_red = np.ones([nstn_prcp, 2])

            for tt in range(ntimes):
            # for tt in range(time1-1, time2): # to accelerate speed
                prea_tar = readata_stn[rr, gg, tt]
                if np.isnan(prea_tar):
                    continue
                if stndata[gg, tt]>0:
                    pstn_tar = 1
                else:
                    pstn_tar = 0
                prea_near = readata_stn[rr, nearloc, tt]
                pstn_near = stndata[nearloc, tt]
                pstn_near[pstn_near > 0] = 1

                # logistic regression
                if np.all(pstn_near == 1):
                    reapop_stn[rr, gg, tt] = 1
                elif np.all(pstn_near == 0) or np.all(prea_near < 0.01):
                    reapop_stn[rr, gg, tt] = 0
                else:
                    x_red[:, 1] = prea_near
                    tx_red = np.transpose(x_red)
                    twx_red = np.matmul(tx_red, w_pcp_red)
                    b = reg.logistic_regression(x_red, twx_red, pstn_near)
                    if np.all(b==0) or np.any(np.isnan(b)):
                        reapop_stn[rr, gg, tt] = np.dot(nearweight, pstn_near)
                    else:
                        zb = - np.dot(np.array([1,prea_tar]), b)
                        reapop_stn[rr, gg, tt] = 1 / (1 + np.exp(zb))
    np.savez_compressed(file_reapop_stn, reapop_stn=reapop_stn, stnlle=stnlle)

########################################################################################################################

if os.path.isfile(file_popmerge_stn):
    datatemp = np.load(file_popmerge_stn)
    bmaweight_stn = datatemp['bmaweight_stn']
    mergepop_stn = datatemp['mergepop_stn']
    bmaweight_grid = datatemp['bmaweight_grid']
    del datatemp
else:
    # estimate bma merging weights for pop
    print('estimate bma merging weights for pop')
    bmaweight_stn = np.nan * np.zeros([12, nstn, reanum], dtype=np.float32)
    for m in range(12):
        print('month',m)
        indm = date_number['mm'] == (m+1)
        for i in range(nstn):
            if np.isnan(stndata[i, 0]):
                continue
            rea = reapop_stn[:, i, indm].T
            obs = stndata[i, indm].copy()
            obs[obs > 0] = 1
            weight, sigma, sigma_s = bma(rea, obs)
            bmaweight_stn[m, i, :] = weight

    # use cross validation to calculate independent merged pop which can be evaluated using station data
    mergepop_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)
    for i in range(nstn):
        if np.mod(i,5000)==0:
            print(i)
        if np.isnan(stndata[i, 0]):
            continue
        nearloc = near_loc_stn[i, :]
        nearweight = near_weight_stn[i, :]
        nearweight = nearweight[nearloc > -1]
        nearweight = nearweight / np.sum(nearweight)
        nearloc = nearloc[nearloc > -1]
        nearweight = np.tile(nearweight,[reanum, 1]).T

        # get bma weight from nearby stations
        weight_i = np.zeros([12, reanum])
        for m in range(12):
            weight_im_near = bmaweight_stn[m, nearloc, :]
            weight_i[m, :] = np.sum(weight_im_near * nearweight, axis=0)

        # merging at the target station
        reapop_merge_i = np.zeros(ntimes)
        for m in range(12):
            indm = date_number['mm'] == (m + 1)
            reapop_stn_im = reapop_stn[:, i, indm]
            weight_im = np.tile(weight_i[m, :], [np.sum(indm), 1]).T
            weight_im[np.isnan(reapop_stn_im)] = np.nan
            reapop_merge_i[indm] = np.nansum(reapop_stn_im * weight_im, axis=0) / np.nansum(weight_im, axis=0)

        mergepop_stn[i, :] = reapop_merge_i

    # interpolate weights to grids
    bmaweight_grid  = np.nan * np.zeros([12, reanum, nrows, ncols], dtype=np.float32)
    for m in range(12):
        for rr in range(reanum):
            bmaweight_grid[m, rr, :, :] = extrapolation(bmaweight_stn[m, :, rr], near_loc_grid, near_dist_grid)
    np.savez_compressed(file_popmerge_stn, bmaweight_stn=bmaweight_stn, bmaweight_grid=bmaweight_grid, mergepop_stn=mergepop_stn)

########################################################################################################################

# this method is quick in speed, but slightly worse than method-1

# # method-2: estimate pop during the proess of QM correction
# date_list, date_number = m_DateList(1979, 2018, 'ByYear')
# rea_pop2 = np.zeros([reanum, nstn, ntimes], dtype=np.float32)
# for rr in range(reanum):
#     for gg in range(nstn):
#         if np.mod(gg,100)==0:
#             print(gg)
#         if np.isnan(stndata[gg, 0]):
#             continue
#         nearloc = near_loc_stn[gg, :]
#         neardist = near_dist_stn[gg, :]
#         nearweight = near_weight_stn[gg, :]
#         neardist = neardist[nearloc > -1]
#         nearweight = nearweight[nearloc > -1]
#         nearweight = nearweight / np.sum(nearweight)
#         nearloc = nearloc[nearloc > -1]
#         nstn_prcp = len(nearloc)
#
#         for mm in range(12):
#             indm = date_number['mm'] == (mm+1)
#             prea_tar = readata_stn[rr, gg, indm]
#             pstn_near = stndata[nearloc,:][:, indm]
#             ecdf_reatar = empirical_cdf(prea_tar, ecdf_prob)
#             popmm_near = np.zeros([nstn_prcp, np.sum(indm)])
#             popmm = np.zeros(np.sum(indm))
#             for i in range(nstn_prcp):
#                 # a simple Qunatile mapping by sorting
#                 ecdf_neari = empirical_cdf(pstn_near[i,:], ecdf_prob)
#                 pqm = cdf_correction(ecdf_prob, ecdf_neari, ecdf_prob, ecdf_reatar, pstn_near[i,:])
#                 pqm[pqm > 0] = 1
#                 popmm_near[i, :] = pqm
#                 popmm = popmm + pqm * nearweight[i]
#             rea_pop2[rr, gg, indm] = popmm

########################################################################################################################

# gridded pop estimation, bma-based merging and error estimation
# year=[1979,2018]
# for y in range(year[0], year[1] + 1):
#     for m in range(12):
for y in range(yearin, yearin + 1):
    for m in range(monthin-1, monthin):
        print('year,month',y,m+1)
        file_reapop = path_pop + '/rea_pop_' + str(y * 100 + m + 1) + '.npz'
        file_bmapop = path_pop + '/bmamerge_pop_' + str(y * 100 + m + 1) + '.npz'
        if os.path.isfile(file_bmapop):
            print('file exists ... continue')
            continue

        # date processing
        indmy = (date_number['yyyy'] == y) & (date_number['mm'] == m + 1)
        mmdays = np.sum(indmy)

        # read raw gridded reanalysis data
        readata_raw = np.nan * np.zeros([reanum, nrows, ncols, mmdays], dtype=np.float32)
        for rr in range(reanum):
            if not (prefix[rr] == 'MERRA2_' and y == 1979):
                filer = path_readowngrid[rr] + '/' + prefix[rr] + 'ds_prcp_' + str(y*100 +m+1) + '.npz'
                d = np.load(filer)
                readata_raw[rr, :, :, :] = d['data']
                del d

        readata_stnym = readata_stn[:, :, indmy].copy()
        stndataym = stndata[:, indmy].copy()

        ################################################################################################################
        print('estimate pop for all grids')
        reapop_grid = np.nan * np.zeros([reanum, nrows, ncols, mmdays], dtype=np.float32)
        if os.path.isfile(file_reapop):
            datatemp = np.load(file_reapop)
            reapop_grid = datatemp['reapop_grid']
            del datatemp
        else:
            for r in range(nrows):
                if np.mod(r, 10) == 0:
                    print(r, nrows)
                for c in range(1000, 1001):
                    if np.isnan(mask[r, c]):
                        continue
                    nearloc = near_loc_grid[r, c, :]
                    neardist = near_dist_grid[r, c, :]
                    nearweight = near_weight_grid[r, c, :]
                    neardist = neardist[nearloc > -1]
                    nearweight = nearweight[nearloc > -1]
                    nearweight = nearweight / np.sum(nearweight)
                    nearloc = nearloc[nearloc > -1]

                    nstn_prcp = len(nearloc)
                    w_pcp_red = np.zeros([nstn_prcp, nstn_prcp])
                    for i in range(nstn_prcp):
                        w_pcp_red[i, i] = nearweight[i]  # eye matrix: stn weight in one-one lien

                    x_red = np.ones([nstn_prcp, 2])
                    for rr in range(reanum):
                        for tt in range(mmdays):
                            prea_tar = readata_raw[rr, r, c, tt]
                            if np.isnan(prea_tar):
                                continue
                            prea_near = readata_stnym[rr, nearloc, tt]
                            pstn_near = stndataym[nearloc, tt]
                            pstn_near[pstn_near > 0] = 1

                            # logistic regression
                            if np.all(pstn_near == 1):
                                reapop_grid[rr, r, c, tt] = 1
                            elif np.all(pstn_near == 0) or np.all(prea_near < 0.01):
                                reapop_grid[rr, r, c, tt] = 0
                            else:
                                x_red[:, 1] = prea_near
                                tx_red = np.transpose(x_red)
                                twx_red = np.matmul(tx_red, w_pcp_red)
                                b = reg.logistic_regression(x_red, twx_red, pstn_near)
                                if np.all(b == 0) or np.any(np.isnan(b)):
                                    reapop_grid[rr, r, c, tt] = np.dot(nearweight, pstn_near)
                                else:
                                    zb = - np.dot(np.array([1, prea_tar]), b)
                                    reapop_grid[rr, r, c, tt] = 1 / (1 + np.exp(zb))

                                # # another choice for logistic regression
                                # # model=SGDClassifier(loss='log')
                                # model = LogisticRegression(solver='lbfgs')
                                # model.fit(np.reshape(prea_near, [-1, 1]), pstn_near, sample_weight=nearweight)
                                # reapop_grid[rr, r, c, tt] = model.predict_proba(np.reshape(prea_tar, [-1, 1]))[0][1]

            np.savez_compressed(file_reapop, reapop_grid=reapop_grid, prefix=prefix)

        ################################################################################################################
        print('Reanalysis merging')
        # start BMA-based merging
        if not os.path.isfile(file_bmapop):
            # initialization
            bma_data = np.nan * np.zeros([nrows, ncols, mmdays], dtype=np.float32)

            # (1) estimate the error of corrected data by interpolating stations
            obs = stndata[:, indmy].copy()
            obs[obs>0]=1
            bma_error = extrapolation(mergepop_stn[:, indmy] - obs, near_loc_grid, near_dist_grid)

            # (2) estimate the value of merged data
            reamerge_weight_gridm = bmaweight_grid[m, :, :, :].copy()
            for i in range(mmdays):
                datai = reapop_grid[:, :, :, i]
                weighti = reamerge_weight_gridm.copy()
                weighti[np.isnan(datai)] = np.nan
                bma_data[:, :, i] = np.nansum(weighti * datai, axis=0) / np.nansum(weighti, axis=0)
            np.savez_compressed(file_bmapop, bma_data=bma_data, bma_error=bma_error,
                                reaname=prefix, latitude=lattar, longitude=lontar)
            del bma_error, bma_data
