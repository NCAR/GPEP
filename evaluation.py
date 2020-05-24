# evaluate downscaled/corrected/merged/regression estimates

import numpy as np
from scipy import io
import auxiliary as au
import os

filestn = '/home/gut428/stndata_whole.npz'
datatemp = np.load(filestn)
prcp_stn = datatemp['prcp_stn']
tmean_stn = datatemp['tmean_stn']
trange_stn = datatemp['trange_stn']
stn_ID=datatemp['stn_ID']
stn_lle=datatemp['stn_lle']
del datatemp
nstn, ntimes = np.shape(prcp_stn)

# 1. evaluate downscaled prcp
file_readown = ['/datastore/GLOBALWATER/CommonData/EMDNA/ERA5_day_ds/ERA5_downto_stn_nearest.npz',
                '/datastore/GLOBALWATER/CommonData/EMDNA/MERRA2_day_ds/MERRA2_downto_stn_nearest.npz',
                '/datastore/GLOBALWATER/CommonData/EMDNA/JRA55_day_ds/JRA55_downto_stn_nearest.npz',
                '/datastore/GLOBALWATER/CommonData/EMDNA/ERA5_day_ds/ERA5_downto_stn_linear.npz',
                '/datastore/GLOBALWATER/CommonData/EMDNA/MERRA2_day_ds/MERRA2_downto_stn_linear.npz',
                '/datastore/GLOBALWATER/CommonData/EMDNA/JRA55_day_ds/JRA55_downto_stn_linear.npz',
                '/datastore/GLOBALWATER/CommonData/EMDNA/ERA5_day_ds/ERA5_downto_stn_GWR.npz',
                '/datastore/GLOBALWATER/CommonData/EMDNA/MERRA2_day_ds/MERRA2_downto_stn_GWR.npz',
                '/datastore/GLOBALWATER/CommonData/EMDNA/JRA55_day_ds/JRA55_downto_stn_GWR.npz']
fileout = 'evaluation_readown.npz'
fileout2 = 'evaluation_readown.mat'
if not os.path.isfile(fileout):
    print('Evaluate downscaled reanalysis')
    fnum = len(file_readown)
    met_prcp = np.zeros([fnum, nstn, 4])
    met_tmean = np.zeros([fnum, nstn, 4])
    met_trange = np.zeros([fnum, nstn, 4])
    for i in range(fnum):
        print('File', i, fnum)
        datatemp = np.load(file_readown[i])
        prcp_readown = datatemp['prcp_readown']
        tmean_readown = datatemp['tmean_readown']
        trange_readown = datatemp['trange_readown']
        if np.shape(prcp_readown)[1]<14610:
            st = 365
        else:
            st = 0

        for g in range(nstn):
            met_prcp[i, g, :] = au.metric(prcp_stn[g, st:], prcp_readown[g, :])
            met_tmean[i, g, :] = au.metric(tmean_stn[g, st:], tmean_readown[g, :])
            met_trange[i, g, :] = au.metric(trange_stn[g, st:], trange_readown[g, :])
        del datatemp

    np.savez_compressed(fileout,met_prcp=met_prcp,met_tmean=met_tmean,met_trange=met_trange,
                        stn_ID=stn_ID,stn_lle=stn_lle,file_readown=file_readown)
    io.savemat(fileout2,{'met_prcp':met_prcp, 'met_tmean':met_tmean, 'met_trange':met_trange,
                         'stn_ID':stn_ID, 'stn_lle':stn_lle,'file_readown':file_readown})