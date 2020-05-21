# main_daily.py produces error estimation at station points for each month
# this script integrate them into one single file and recover regression estimates from error estimate

import numpy as np
import calendar
from auxiliary_merge import *
import auxiliary as au

path = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout'
outfile = '/home/gut428/daily_regression_stn.npz'
nstn = 27275
ntimes = 14610
date_list, date_number = m_DateList(1979, 2018, 'ByYear')

prcp = np.zeros([nstn, ntimes], dtype=np.float32)
prcp_bc = np.zeros([nstn, ntimes], dtype=np.float32)
tmean = np.zeros([nstn, ntimes], dtype=np.float32)
trange = np.zeros([nstn, ntimes], dtype=np.float32)

for y in range(1979, 2019):
    for m in range(1, 13):
        print('Integrate regression data: Year month:', y, m)
        indym = (date_number['yyyy'] == y) & (date_number['mm'] == m)

        date_cal_start = y * 10000 + m * 100 + 1
        date_cal_end = y * 10000 + m * 100 + calendar.monthrange(y, m)[1]
        datestr = str(date_cal_start) + '-' + str(date_cal_end)
        file1 = path + '/error_' + datestr + '.npz'
        file2 = path + '/stndata_' + datestr + '.npz'
        d1 = np.load(file1)
        d2 = np.load(file2)

        temp = d2['prcp_stn_daily']
        prcp[:, indym] = temp + d1['prcp_err_stn_raw']
        prcp_bc[:, indym] = au.retransform(au.transform(temp, 4, 'box-cox') + d1['prcp_err_stn_bc'], 4, 'box-cox')
        tmean[:, indym] = d2['tmean_stn_daily'] + d1['tmean_err_stn']
        trange[:, indym] = d2['trange_stn_daily'] + d1['trange_err_stn']

np.savez_compressed(outfile, prcp=prcp, prcp_bc=prcp_bc, tmean=tmean, trange=trange)
