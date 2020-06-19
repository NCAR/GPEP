# check the quality of station data based on regression estimates to ensure that some toxic stations are excluded
# note: s2_ and s3_ must be re-run after executing this step because stations are changed

import numpy as np

nearnum = 10  # compare target station with its nearby stations to determine climatological anomaly

filestn = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_original.npz'
filereg = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_original/regression_stn.npz'
file_nearstn = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_original/nearstn_catalog.npz'
# filestn = '/Users/localuser/Research/EMDNA/stndata_whole.npz'
# filereg = '/Users/localuser/Research/EMDNA/regression/daily_regression_stn.npz'
# file_nearstn = '/Users/localuser/Research/EMDNA/regression/weight_nearstn.npz'

# outfiles
filestn_new = '/home/gut428/stndata_aftercheck.npz'

########################################################################################################################

# precipitation
var = 'prcp'
d = np.load(filestn)
dstn = d[var + '_stn']
d = np.load(filereg)
dreg = d[var]

d = np.load(file_nearstn)
if var == 'prcp':
    nearloc = d['near_stn_prcpLoc']
else:
    nearloc = d['near_stn_tempLoc']

# calculate two types of ratio
nstn, ntimes = np.shape(dstn)
dstnm = np.nanmean(dstn, axis=1)
dregm = np.nanmean(dreg, axis=1)
ratio_vsreg = np.abs(dstnm) / np.abs(dregm)
ratio_vsnear = np.nan * np.zeros(nstn)
for i in range(nstn):
    if np.isnan(dstnm[i]):
        continue
    zz = nearloc[i, 0:nearnum]
    zz = zz[zz > -1]
    ratio_vsnear[i] = abs(dstnm[i] / np.nanmean(dstnm[zz]))

# determine toxic stations: two-tailed percentile check
# the thresholds are determined by trial-and-error
# for precipitation, stations with too dry precipitation is more common. besides, we prefer to keep more stations
# with heavy precipitation to benefit extreme precipitation estimation
p1 = np.nanpercentile(ratio_vsreg, 99.9)
p2 = np.nanpercentile(ratio_vsnear, 99.9)
p3 = np.nanpercentile(ratio_vsreg, 1)
p4 = np.nanpercentile(ratio_vsnear, 1)

# toxic_prcp=np.where(((ratio_vsreg<p3) & (ratio_vsnear<p4)))[0]
toxic_prcp = np.where(((ratio_vsreg > p1) & (ratio_vsnear > p2)) |
                      ((ratio_vsreg < p3) & (ratio_vsnear < p4)))[0]

del dstn, dregm

########################################################################################################################

# tmean/trange
vart = ['tmean', 'trange']
toxic_temp = [' '] * 2
for v in range(len(vart)):
    var = vart[v]
    d = np.load(filestn)
    dstn = d[var + '_stn']
    d = np.load(filereg)
    dreg = d[var]

    d = np.load(file_nearstn)
    if var == 'prcp':
        nearloc = d['near_stn_prcpLoc']
    else:
        nearloc = d['near_stn_tempLoc']

    nstn, ntimes = np.shape(dstn)
    dstnm = np.nanmean(dstn, axis=1)
    dregm = np.nanmean(dreg, axis=1)
    diff_vsreg = np.abs(dstnm - dregm)
    diff_vsnear = np.nan * np.zeros(nstn)
    for i in range(nstn):
        if np.isnan(dstnm[i]):
            continue
        zz = nearloc[i, 0:nearnum]
        zz = zz[zz > -1]
        diff_vsnear[i] = np.abs(dstnm[i] - np.mean(dstnm[zz]))

    p1 = np.nanpercentile(diff_vsreg, 99)
    p2 = np.nanpercentile(diff_vsnear, 95)
    toxic_temp[v] = np.where((diff_vsreg > p1) & (diff_vsnear > p2))[0]
toxic_temp = np.union1d(toxic_temp[0], toxic_temp[1])

########################################################################################################################

# save station data in a new file
d = np.load(filestn)
prcp_stn = d['prcp_stn']
tmean_stn = d['tmean_stn']
trange_stn = d['trange_stn']
date_ymd = d['date_ymd']
stnID = d['stnID']
stninfo = d['stninfo']

prcp_stn[toxic_prcp, :] = np.nan
tmean_stn[toxic_temp, :] = np.nan
trange_stn[toxic_temp, :] = np.nan

np.savez_compressed(filestn_new, prcp_stn=prcp_stn, tmean_stn=tmean_stn, trange_stn=trange_stn,
                    date_ymd=date_ymd, stnID=stnID, stninfo=stninfo, toxic_prcp=toxic_prcp, toxic_temp=toxic_temp)
