import numpy as np
from auxiliary_merge import m_DateList
vars=['prcp','tmean','trange']
reanum=3
nstn=27275
ntimes=14610
date_list, date_number = m_DateList(1979, 2018, 'ByYear')

for i in range(len(vars)):
    file_corrmerge_stn='/home/gut428/ReanalysisCorrMerge/Reanalysis_merge/mergecorr_stn_' + vars[i] + '_GWRQM_BMA.npz'

    reamerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)
    reamerge_weight_stn = np.nan * np.zeros([12, nstn, reanum], dtype=np.float32)
    reacorr_stn = np.nan * np.zeros([reanum, nstn, ntimes], dtype=np.float32)
    # for each month
    for m in range(12):
        print('month', m + 1)
        filem = '/home/gut428/ReanalysisCorrMerge/Reanalysis_merge/mergecorr_stn_' \
                + vars[i] + '_GWRQM_BMA_' + str(m+1) + '.npz'
        indm = date_number['mm'] == (m + 1)

        d=np.load(filem)
        date_ymd=d['date_ymd']
        prefix=d['prefix']
        stninfo=d['stninfo']
        reamerge_stn[:, indm] = d['reamerge_stn'][:, indm]
        reacorr_stn[:, :, indm] = d['reacorr_stn'][:, :, indm]
        reamerge_weight_stn[m, :, :] = d['reamerge_weight_stn'][m, :, :]

    # the variables are independent with their concurrent stations. thus, station data can be used to evaluate them
    np.savez_compressed(file_corrmerge_stn, reamerge_stn=reamerge_stn, reamerge_weight_stn=reamerge_weight_stn,
                        reacorr_stn=reacorr_stn, date_ymd=date_ymd, prefix=prefix, stninfo=stninfo)