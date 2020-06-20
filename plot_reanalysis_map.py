import numpy as np
from matplotlib import pyplot as plt
import h5py
import datetime

year = [1979, 2018]
month = [1, 12]
nrows = 800
ncols = 1300
vars=['prcp','tmean','trange']
clims=[[0,10], [-40,40], [0,20]]
reaname=['ERA5','MERRA2','JRA55']

for v in range(len(vars)):
    var = vars[v]
    for rea in reaname:
        varname = 'data'
        outpath = '/home/gut428/figures/' + rea + '_day_ds'
        inpath = '/datastore/GLOBALWATER/CommonData/EMDNA_new/' + rea + '_day_ds'
        # monthly
        for y in range(2018, 2019):
            print('year', y)
            for m in range(month[0], month[1] + 1):
                file = inpath + '/' + rea + '_ds_' + var + '_' + str(y * 100 + m) + '.npz'
                figm = outpath + '/' + var + '_' + str(y * 100 + m) + '.png'
                d = np.load(file)
                data = d[varname]
                data = np.mean(data, axis=2)
                plt.imshow(data, cmap='nipy_spectral')
                plt.colorbar(orientation='horizontal', shrink=0.6, pad=0.05)
                plt.clim(clims[v])
                plt.title('Date:'+str(y * 100 + m))
                plt.axis('off')
                plt.savefig(figm, dpi=600)
                plt.close()

        # monthly climatology
        for m in range(month[0], month[1] + 1):
            print('month',m)
            nday = 0
            data = np.zeros([nrows, ncols])
            figm = outpath + '/' + var + '_' + str(m) + '.png'
            for y in range(year[0], year[1] + 1):
                file = inpath + '/' + rea + '_ds_' + var + '_' + str(y * 100 + m) + '.npz'
                d = np.load(file)
                dy = d[varname]
                data = data + np.sum(dy, axis=2)
                nday = nday + np.shape(dy)[2]
            data = data / nday
            plt.imshow(data, cmap='nipy_spectral')
            plt.colorbar(orientation='horizontal', shrink=0.6, pad=0.05)
            plt.clim(clims[v])
            plt.axis('off')
            plt.savefig(figm, dpi=600)
            plt.close()


# raw reanalysis data
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

for y in range(2018,2019):
    date_list, date_number = m_DateList(y, y, 'ByYear')
    for rea in reaname:
        # prcp
        infile = '/datastore/GLOBALWATER/CommonData/EMDNA/' + rea + '_day_raw/' + rea + '_prcp_' + str(y) + '.mat'
        outpath = '/home/gut428/figures/' + rea + '_day_raw'
        datatemp = {}
        f = h5py.File(infile, 'r')
        for k, v in f.items():
            datatemp[k] = np.array(v)
        dataori = datatemp['data']
        dataori = np.transpose(dataori, [2, 1, 0])
        del datatemp
        f.close()
        for m in range(12):
            figm = outpath + '/prcp_' + str(y*100+m+1) + '.png'
            indm = date_number['mm'] == m + 1
            data = np.mean(dataori[:,:,indm], axis=2)
            plt.imshow(data, cmap='nipy_spectral')
            plt.colorbar(orientation='horizontal', shrink=0.6, pad=0.05)
            plt.clim(clims[0])
            plt.title('Date:' + str(y * 100 + m+1))
            plt.axis('off')
            plt.savefig(figm, dpi=600)
            plt.close()
        # tmean/trange
        infile = '/datastore/GLOBALWATER/CommonData/EMDNA/' + rea + '_day_raw/' + rea + '_tmin_' + str(y) + '.mat'
        datatemp = {}
        f = h5py.File(infile, 'r')
        for k, v in f.items():
            datatemp[k] = np.array(v)
        tmin = datatemp['data']
        tmin = np.transpose(tmin, [2, 1, 0])
        del datatemp
        f.close()
        infile = '/datastore/GLOBALWATER/CommonData/EMDNA/' + rea + '_day_raw/' + rea + '_tmax_' + str(y) + '.mat'
        datatemp = {}
        f = h5py.File(infile, 'r')
        for k, v in f.items():
            datatemp[k] = np.array(v)
        tmax = datatemp['data']
        tmax = np.transpose(tmax, [2, 1, 0])
        del datatemp
        f.close()
        tmean=(tmin+tmax)/2
        trange=np.abs(tmax-tmin)
        for m in range(12):
            figm = outpath + '/tmean_' + str(y*100+m+1) + '.png'
            indm = date_number['mm'] == m + 1
            data = np.mean(tmean[:, :, indm], axis=2)
            plt.imshow(data, cmap='nipy_spectral')
            plt.colorbar(orientation='horizontal', shrink=0.6, pad=0.05)
            plt.clim(clims[1])
            plt.title('Date:' + str(y * 100 + m+1))
            plt.axis('off')
            plt.savefig(figm, dpi=600)
            plt.close()
            figm = outpath + '/trange_' + str(y*100+m+1) + '.png'
            indm = date_number['mm'] == m + 1
            data = np.mean(trange[:, :, indm], axis=2)
            plt.imshow(data, cmap='nipy_spectral')
            plt.colorbar(orientation='horizontal', shrink=0.6, pad=0.05)
            plt.clim(clims[2])
            plt.title('Date:' + str(y * 100 + m+1))
            plt.axis('off')
            plt.savefig(figm, dpi=600)
            plt.close()