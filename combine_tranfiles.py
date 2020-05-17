import numpy as np
import calendar
import os

path_trans = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout'
path_notrans = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout_notrans'
path_out = '/home/gut428/PyGMETout'
year = [1979,2018]

for y in range(year[0],year[1]+1):
    for m in range(12):
        date_cal_start = y * 10000 + (m + 1) * 100 + 1
        date_cal_end = y * 10000 + (m + 1) * 100 + calendar.monthrange(y, m + 1)[1]
        print('Date', date_cal_start, date_cal_end)
        datestr = str(date_cal_start) + '-' + str(date_cal_end)

        # combine two outout files
        fileo1 = path_trans + '/output_' + datestr + '.npz'
        fileo2 = path_notrans + '/output_notrans_' + datestr + '.npz'
        fileoout = path_out + '/output_' + datestr + '.npz'
        fileoout2 = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout_tran_notran/output_' + datestr + '.npz'
        if (not os.path.isfile(fileoout)) and (not os.path.isfile(fileoout2)):
            d1 = np.load(fileo1)
            d2 = np.load(fileo2)
            pop = d1['pop']
            pcp_bc = d1['pcp']
            pcp_err_bc = d1['pcp_err']
            y_max_bc = d1['y_max']
            tmean = d1['tmean']
            trange = d1['trange']
            tmean_err = d1['tmean_err']
            trange_err = d1['trange_err']
            mean_autocorr_daily = d1['mean_autocorr_daily']
            mean_tp_corr_daily = d1['mean_tp_corr_daily']
            pcp_raw = d2['pcp']
            pcp_err_raw = d2['pcp_err']
            y_max_raw = d2['y_max']

            np.savez_compressed(fileoout, pop=pop, pcp_bc=pcp_bc, pcp_err_bc=pcp_err_bc, y_max_bc=y_max_bc,
                                pcp_raw=pcp_raw, pcp_err_raw=pcp_err_raw, y_max_raw=y_max_raw,
                                tmean=tmean, trange=trange, tmean_err=tmean_err, trange_err=trange_err,
                                mean_autocorr_daily=mean_autocorr_daily, mean_tp_corr_daily=mean_tp_corr_daily)

        # combine two error files
        filee1 = path_trans + '/error_' + datestr + '.npz'
        filee2 = path_notrans + '/error_notrans_' + datestr + '.npz'
        fileeout = path_out + '/error_' + datestr + '.npz'
        fileeout2 = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout_tran_notran/error_' + datestr + '.npz'
        if (not os.path.isfile(fileeout)) and (not os.path.isfile(fileeout2)):
            d1 = np.load(filee1)
            d2 = np.load(filee2)
            pcp_err_stn_bc = d1['pcp_err_stn']
            tmean_err_stn = d1['tmean_err_stn']
            trange_err_stn = d1['trange_err_stn']
            stninfo = d1['stninfo']
            pcp_err_stn_raw = d2['pcp_err_stn']

            np.savez_compressed(fileeout, pcp_err_stn_bc=pcp_err_stn_bc, pcp_err_stn_raw=pcp_err_stn_raw,
                                tmean_err_stn=tmean_err_stn, trange_err_stn=trange_err_stn, stninfo=stninfo)