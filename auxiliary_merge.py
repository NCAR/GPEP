import datetime
import numpy as np
import sys
import auxiliary as au

def m_DateList(year_start, year_end, mode):
    # generate a date list (yyyymmdd) between start year and end year
    # mode: 'ByDay', 'ByMonth', 'ByYear': time scales of input files
    date_start = datetime.date(year_start, 1, 1)
    date_end = datetime.date(year_end, 12, 31)
    daynum = (date_end - date_start).days + 1

    # generate date in format: yyyymmdd
    date_ymd = np.zeros(daynum, dtype=int)
    dated = date_start
    for d in range(daynum):
        if d > 0:
            dated = dated + datetime.timedelta(days=1)
        date_ymd[d] = int(dated.strftime("%Y%m%d"))
    date_number = {'yyyymmdd': date_ymd,
                   'yyyymm': np.floor(date_ymd / 100).astype(int),
                   'yyyy': np.floor(date_ymd / 10000).astype(int),
                   'mm': np.floor(np.mod(date_ymd, 10000) / 100).astype(int),
                   'dd': np.mod(date_ymd, 100).astype(int)}

    # generate file list
    if mode == 'ByDay':
        datemode = date_number['yyyymmdd']
    else:
        if mode == 'ByMonth':
            datemode = date_number['yyyymm']
        elif mode == 'ByYear':
            datemode = date_number['yyyy']
        datemode = np.unique(datemode)

    date_list = [' '] * len(datemode)
    for i in range(len(datemode)):
        date_list[i] = str(datemode[i])

    return date_list, date_number

def box_cox_transform(data, exp=0.25):
    return (data ** exp - 1) / exp

def box_cox_recover(data, exp=0.25):
    dataout = (data * exp + 1) ** (1/exp)
    dataout[data < -1/exp] = 0
    return dataout

def calmetric(dtar, dref, metname='RMSE'):
    if np.ndim(dtar) == 1:
        dtar = dtar[np.newaxis, :]
        dref = dref[np.newaxis, :]
    nstn, ntimes = np.shape(dtar)
    metout = np.nan * np.zeros(nstn, dtype=np.float32)
    if metname == 'RMSE':
        for i in range(nstn):
            metout[i] = np.sqrt(np.nansum(np.square(dtar[i, :] - dref[i, :])) / ntimes)  # RMSE
    elif metname == 'CC':
        for i in range(nstn):
            temp = np.corrcoef(dtar[i, :], dref[i, :])
            metout[i] = temp[0, 1]
    else:
        sys.exit('Unkown metric name')
    return metout


def ismember(a, b):
    # tf = np.in1d(a,b) # for newer versions of numpy
    tf = np.array([i in b for i in a])
    u = np.unique(a[tf])
    index = np.array([(np.where(b == i))[0][-1] if t else 0 for i, t in zip(a, tf)])
    return tf, index

def extrapolation(datain, nearstn_loc, nearstn_dist):
    # datain: one or multiple time steps
    wexp = 3
    if np.ndim(datain) == 1:  # add time axis
        datain = datain[:, np.newaxis]

    if np.ndim(nearstn_loc) == 2:  # extrapolate to station points
        num = np.shape(nearstn_loc)[0]
        ntimes = np.shape(datain)[1]
        dataout = np.nan * np.zeros([num, ntimes], dtype=np.float32)
        for i in range(num):
            if not nearstn_loc[i, 0]>=0:
                continue
            nearloci = nearstn_loc[i, :]
            indloci = nearloci > -1
            dataini = datain[nearloci[indloci], :]
            disti = nearstn_dist[i, indloci]
            weighti = au.distanceweight(disti, np.max(disti) + 1, wexp)
            weighti[np.isnan(dataini[:,0])] = np.nan
            weighti = weighti / np.nansum(weighti)
            weighti2 = np.tile(weighti,[ntimes,1]).T
            dataout[i, :] = np.nansum(dataini * weighti2, axis=0)
    elif np.ndim(nearstn_loc) == 3:  # extrapolate to gridds
        nrows, ncols, nearnum = np.shape(nearstn_loc)
        nstn, ntimes = np.shape(datain)
        dataout = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
        for r in range(nrows):
            for c in range(ncols):
                if not nearstn_loc[r, c, 0]>=0:
                    continue
                nearloci = nearstn_loc[r, c, :]
                indloci = nearloci > -1
                dataini = datain[nearloci[indloci], :]
                disti = nearstn_dist[r, c, indloci]
                weighti = au.distanceweight(disti, np.max(disti) + 1, wexp)
                weighti[np.isnan(dataini[:,0])]=np.nan
                weighti = weighti / np.nansum(weighti)
                weighti2 = np.tile(weighti, [ntimes, 1]).T
                dataout[r, c, :] = np.nansum(dataini * weighti2, axis=0)
    else:
        print('The dimensions of tarlat or tarlon are larger than 2')
        sys.exit()
    if ntimes == 1:
        dataout = np.squeeze(dataout)
    return dataout
