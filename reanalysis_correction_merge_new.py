import numpy as np
import os, sys
import auxiliary as au

def empirical_cdf(data, probtar):
    # data: vector of data
    data2 = data[~np.isnan(data)]
    if len(data2) > 0:
        ds = np.sort(data)
        probreal = np.arange(len(data2)) / (len(data2)+1)
        ecdf_out = np.interp(probtar, probreal, ds)
    else:
        ecdf_out = np.nan * np.zeros(len(probtar))
    return ecdf_out

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


# settings
gmet_stndatafile = '/home/gut428/stndata_whole.npz'
outpath_ecdf_stn = ''
bin = 500
vars = ['prcp','tmean','trange']

########################################################################################################################
# basic information

datatemp = np.load(gmet_stndatafile)
stnlle = datatemp['stn_lle']
nstn = np.shape(stnlle)[0]

########################################################################################################################

# find near stations for all grids and station

if var == 'trange':
    vari = 'tmean'  # trange and tmean have the same index
else:
    vari = var

near_stnfile = near_path + '/near_stn_' + vari + '.npz'
near_gridfile = near_path + '/near_grid_' + vari + '.npz'

if os.path.isfile(near_stnfile):
    print('load near station information for points')
    with np.load(near_stnfile) as datatemp:
        nearstn_locl1 = datatemp['nearstn_locl1']
        nearstn_distl1 = datatemp['nearstn_distl1']
        nearstn_locl2 = datatemp['nearstn_locl2']
        nearstn_distl2 = datatemp['nearstn_distl2']
    del datatemp
else:
    print('find near stations for station points')
    # layer-1
    nstn_testl1 = np.shape(taintestindex[vari + '_testindex1'])[1]
    nearstn_locl1 = -1 * np.ones([dividenum, nstn_testl1, nearnum], dtype=int)
    nearstn_distl1 = -1 * np.ones([dividenum, nstn_testl1, nearnum], dtype=np.float32)
    for lay1 in range(dividenum):
        trainindex1 = taintestindex[vari + '_trainindex1'][lay1, :]
        testindex1 = taintestindex[vari + '_testindex1'][lay1, :]
        nearstn_locl1[lay1, :, :], nearstn_distl1[lay1, :, :] \
            = findnearstn(stnlle[trainindex1, 0], stnlle[trainindex1, 1],
                          stnlle[testindex1, 0], stnlle[testindex1, 1], nearnum, 0)
    # layer-2
    nstn_testl2 = np.shape(taintestindex[vari + '_testindex2'])[2]
    nearstn_locl2 = -1 * np.ones([dividenum, dividenum, nstn_testl2, nearnum], dtype=int)
    nearstn_distl2 = -1 * np.ones([dividenum, dividenum, nstn_testl2, nearnum], dtype=np.float32)
    for lay1 in range(dividenum):
        for lay2 in range(dividenum):
            trainindex2 = taintestindex[vari + '_trainindex2'][lay1, lay2, :]
            testindex2 = taintestindex[vari + '_testindex2'][lay1, lay2, :]
            nearstn_locl2[lay1, lay2, :, :], nearstn_distl2[lay1, lay2, :, :] \
                = findnearstn(stnlle[trainindex2, 0], stnlle[trainindex2, 1],
                              stnlle[testindex2, 0], stnlle[testindex2, 1], nearnum, 0)

    np.savez_compressed(near_stnfile, nearstn_locl1=nearstn_locl1, nearstn_distl1=nearstn_distl1,
                        nearstn_locl2=nearstn_locl2, nearstn_distl2=nearstn_distl2)

if os.path.isfile(near_gridfile):
    print('load near station information for grids')
    with np.load(near_gridfile) as datatemp:
        neargrid_loc = datatemp['neargrid_loc']
        neargrid_dist = datatemp['neargrid_dist']
else:
    print('find near stations for grids')
    stnlle_in = stnlle.copy()
    stnlle_in[np.isnan(stndata[:, 0]), 0:2] = np.nan
    neargrid_loc, neargrid_dist = findnearstn(stnlle_in[:, 0], stnlle_in[:, 1], lattarm, lontarm, nearnum, 0)
    np.savez_compressed(near_gridfile,neargrid_loc=neargrid_loc,neargrid_dist=neargrid_dist)



########################################################################################################################

# obtain empirical CDF curves for all stations

# get probability vector
prob = np.arange(0, 1 + 1 / bin, 1 / bin)

# loop for each var and each station
datatemp = np.load(gmet_stndatafile)
stnlle = datatemp['stn_lle']
nstn = np.shape(stnlle)[0]
for var in vars:
    print('ECDF var:', var)
    filev = outpath_ecdf_stn + '/ecdf_station_' + var + '.npz'

    if os.path.isfile(filev):
        print('filev exists...')
        continue

    ecdf = np.zeros([nstn, bin + 1], dtype=np.float32)
    stndata = datatemp[var + '_stn']
    for i in range(nstn):
        ecdf[i, :] = empirical_cdf(stndata[i, :], prob)

    np.savez_compressed(filev,ecdf=ecdf, prob=prob, stnlle=stnlle)

del datatemp

########################################################################################################################

# correct reanalysis data at each station
vars = ['prcp','tmean','trange']
file_readownstn = ['/datastore/GLOBALWATER/CommonData/EMDNA/ERA5_day_ds/ERA5_downto_stn.npz',
                   '/datastore/GLOBALWATER/CommonData/EMDNA/MERRA2_day_ds/MERRA2_downto_stn.npz',
                   '/datastore/GLOBALWATER/CommonData/EMDNA/JRA55_day_ds/JRA55_downto_stn.npz']

for var in vars:
    print('ECDF correction:', var)

    # load station data
    datatemp = np.load(gmet_stndatafile)
    stndata = datatemp[var + '_stn']
    del datatemp

    # load station ecdf
    filev = outpath_ecdf_stn + '/ecdf_station_' + var + '.npz'
    datatemp = np.load(filev)
    ecdf_stn = datatemp['ecdf']
    prob_stn = datatemp['prob']
    del datatemp

    for filer in file_readownstn:
        # load reanalysis data
        datatemp = np.load(filer)
        readown = datatemp[var + '_readown']
        del datatemp

        # qm correction
        reacorr = np.zeros(np.shape(readown), dtype=np.float32)
        # ntimes = np.shape(reacorr)[1]
        # prob_rea = np.arange(0, 1 + 1 / ntimes, 1 / ntimes)
        prob_rea = prob_stn
        for i in range(nstn):
            ecdf_rea = empirical_cdf(readown[i, :], prob_rea)




reanum = len(file_readownstn)

print('load downscaled reanalysis data at station points')
readata_stn = np.nan * np.zeros([reanum, nstn, ntimes], dtype=np.float32)
for rr in range(reanum):
    dr = np.load(file_readownstn[rr])
    temp = dr[var + '_readown']
    if prefix[rr] == 'MERRA2_':  # unify the time length of all data as MERRA2 lacks 1979
        add = np.nan * np.zeros([nstn, 365])
        temp = np.concatenate((add, temp), axis=1)
    readata_stn[rr, :, :] = temp
    del dr, temp
if var == 'prcp':
    readata_stn[readata_stn<0] = 0