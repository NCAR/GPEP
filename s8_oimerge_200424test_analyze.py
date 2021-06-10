import numpy as np
from scipy import io

def cal_commonmetric(d1, d2):
    # d2 is obs, d1 is estimate
    ind = ~np.isnan(d1) & ~np.isnan(d2)
    d1 = d1[ind]
    d2 = d2[ind]
    cc = np.corrcoef(d1,d2)
    cc = cc[0,1]
    rmse = np.sqrt( np.nanmean((d1-d2)**2) )
    me = np.nanmean(d1-d2)
    kge = KGE_C2020(d2, d1)
    return np.array([cc, rmse, me, kge['KGE'], kge['r'], kge['alpha'], kge['beta']])

def KGE_C2020(Obs, Pre):
    # KGE'' proposed by Clark Martyn et al 2020
    ind_nan = np.isnan(Obs) | np.isnan(Pre)
    Obs = Obs[~ind_nan]
    Pre = Pre[~ind_nan]
    if len(Obs)>2:
        pre_mean = np.mean(Pre, axis=0, dtype=np.float64)
        obs_mean = np.mean(Obs, axis=0, dtype=np.float64)
        r = np.sum((Pre - pre_mean) * (Obs - obs_mean), axis=0, dtype=np.float64) / \
            np.sqrt(np.sum((Pre - pre_mean) ** 2, axis=0, dtype=np.float64) *
                    np.sum((Obs - obs_mean) ** 2, dtype=np.float64))
        # calculate error in spread of flow alpha
        alpha = np.std(Pre, axis=0) / np.std(Obs, dtype=np.float64)
        # calculate error in volume beta (bias of mean discharge)
        beta = (np.mean(Pre) - np.mean(Obs)) / np.std(Obs)
        # calculate the Kling-Gupta Efficiency KGE
        KGE = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta) ** 2)
        KGEgroup = {'KGE': KGE, 'r': r, 'alpha': alpha, 'beta': beta}
    else:
        KGEgroup = {'KGE': np.nan, 'r': np.nan,
                    'alpha': np.nan, 'beta': np.nan}
    return KGEgroup


file_stn = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_aftercheck.npz'
file_oinew = '/home/gut428/OImerge_prcp_loo_20200424/OImerge_stn_GWRBMA_prcp.npz'
file_oiold = '/datastore/GLOBALWATER/CommonData/EMDNA_new/OImerge_GWRLSBMA/OImerge_stn_GWRLSBMA_prcp.npz'


# load data
datatemp = np.load(file_stn)
stndata = datatemp['prcp_stn']
stninfo = datatemp['stninfo']
stnID = datatemp['stnID']
date_ymd = datatemp['date_ymd']
nstn, ntimes = np.shape(stndata)
del datatemp
date_yyyy = (date_ymd/10000).astype(int)
date_mm = (np.mod(date_ymd, 10000)/100).astype(int)


oimerge_old = np.load(file_oiold)
oimerge_old = oimerge_old['oimerge_stn']

oimerge_new = np.load(file_oinew)
oimerge_new = oimerge_new['oimerge_stn']

# evaluate
met_old = np.nan * np.zeros([nstn, 7])
met_new = np.nan * np.zeros([nstn, 7])
for i in range(nstn):
    if np.mod(i, 1000) == 0:
        print(i)
    met_old[i] = cal_commonmetric(oimerge_old[i], stndata[i])
    met_new[i] = cal_commonmetric(oimerge_new[i], stndata[i])

print(np.nanmean(met_old, axis=0))
print(np.nanmean(met_new, axis=0))
print(np.nanmedian(met_old, axis=0))
print(np.nanmedian(met_new, axis=0))
