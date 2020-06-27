# prepare data for GMET probabilistic estimation based on OI merging and original station data
import numpy as np
import auxiliary as au
from calendar import monthrange
import datetime as dt
import netCDF4 as nc
import sys

year=int(sys.argv[1])

# control parameters
yearall = [year, year]
monthall = [1,12]

# ### Mac settings
# path_oi = '/Users/localuser/Research/EMDNA/oimerge'
# near_file_GMET = '/Users/localuser/Research/EMDNA/regression/weight_nearstn.npz'
# gmet_stndatafile = '/Users/localuser/Research/EMDNA/stndata_whole.npz'
# FileStnInfo = '/Users/localuser/Research/EMDNA/basicinfo/stnlist_whole.txt'  # station basic information (lists)
# FileGridInfo = '/Users/localuser/Research/EMDNA/basicinfo/gridinfo_whole.nc'  # study area information
# ### Mac settings

### Plato settings
path_oi = '/datastore/GLOBALWATER/CommonData/EMDNA_new/OImerge_GWRLSBMA'
path_oi_pop = '/datastore/GLOBALWATER/CommonData/EMDNA_new/OImerge_pop'
near_file_GMET = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck/nearstn_catalog.npz'
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_aftercheck.npz'
FileStnInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/stnlist_whole.txt'  # station basic information (lists)
FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/gridinfo_whole.nc'  # study area information
### Plato settings
outpath = '/home/gut428/GMET_OIinput'

########################################################################################################################

print('Read study area basic information')
# load near station information
datatemp = np.load(near_file_GMET)
near_grid_prcpLoc = datatemp['near_grid_prcpLoc']
# near_grid_tempLoc = datatemp['near_grid_tempLoc']
near_grid_prcpLoc = np.flipud(near_grid_prcpLoc)
# near_grid_tempLoc = np.flipud(near_grid_tempLoc)

# station location and attribute information
# stninfo: [ stations, 1/lat/lon/elev/slope_ns/slope_we ]
stnID, stninfo = au.readstnlist(FileStnInfo)

# meshed lat/lon of the target region
ncfid = nc.Dataset(FileGridInfo)
lattarm = ncfid.variables['latitude'][:].data
lattarm = np.flipud(lattarm)
lontarm = ncfid.variables['longitude'][:].data
ncfid.close()
lontar = lontarm[0, :]
lattar = lattarm[:, 0]
nrows = len(lattar)
ncols = len(lontar)

# grid infomation
ncfid = nc.Dataset(FileGridInfo)
gridlat = ncfid.variables['latitude'][:].data
gridlon = ncfid.variables['longitude'][:].data
gridele = ncfid.variables['elev'][:].data
gridgns = ncfid.variables['gradient_n_s'][:].data
gridgwe = ncfid.variables['gradient_w_e'][:].data
mask = ncfid.variables['mask'][:].data  # 1: grids to be considered; the other values: invalid grids
ncfid.close()

gridinfo = np.zeros([nrows, ncols, 6])
gridinfo[:, :, 0] = 1
gridinfo[:, :, 1] = gridlat
gridinfo[:, :, 2] = gridlon
gridinfo[:, :, 3] = gridele
gridinfo[:, :, 4] = gridgns
gridinfo[:, :, 5] = gridgwe
del gridlat, gridlon, gridele, gridgns, gridgwe

# date information
datatemp = np.load(gmet_stndatafile)
prcpstn = datatemp['prcp_stn']
tmeanstn = datatemp['tmean_stn']
trangestn = datatemp['trange_stn']
date_ymd = datatemp['date_ymd']
del datatemp
date_yyyy = (date_ymd/10000).astype(int)
date_mm = (np.mod(date_ymd, 10000)/100).astype(int)

########################################################################################################################

for year in range(yearall[0],yearall[1]+1):
    for month in range(monthall[0], monthall[1]+1):
        print(year, month)
        # output files
        fileout = outpath  + '/reg_' + str(year*100+month) + '.nc'
        # if os.path.isfile(fileout):
        #     print('file exists')
        #     continue
        indym = (date_yyyy == year) & (date_mm == month)

        ################################################################################################################

        # time
        date_cal_start = year * 10000 + month * 100 + 1
        date_cal_end = year * 10000 + month * 100 + monthrange(year, month)[1]
        date_cal_start2 = dt.datetime.strptime(str(date_cal_start), '%Y%m%d')
        date_cal_end2 = dt.datetime.strptime(str(date_cal_end), '%Y%m%d')
        ntimes = (date_cal_end2 - date_cal_start2).days + 1  # time steps to be processed

        # seconds since 1970-1-1 0:0:0
        daydiff = (date_cal_start2 - dt.datetime(1970, 1, 1)).days
        seconds = (np.arange(ntimes) + daydiff) * 86400

        # datelist: yyyymmdd
        yyyymmdd = np.zeros(ntimes, dtype=int)
        for d in range(ntimes):
            dated = date_cal_start2 + dt.timedelta(days=d)
            yyyymmdd[d] = int(dated.strftime("%Y%m%d"))
        yyyymm = np.floor(yyyymmdd / 100).astype(int)
        mm = np.floor(np.mod(yyyymmdd, 10000) / 100).astype(int)

        ################################################################################################################

        # derive max value for each grid pixels
        # load station data
        prcpstn_ym = prcpstn[:, indym]
        # tmeanstn_ym = prcpstn[:, indym]
        # trangestn_ym = prcpstn[:, indym]

        prcp_max = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
        # tmean_max = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
        # trange_max = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
        # prcp_min = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
        # tmean_min = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
        # trange_min = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
        for r in range(nrows):
            for c in range(ncols):
                if near_grid_prcpLoc[r, c, 0] > -1:
                    nearloci = near_grid_prcpLoc[r, c, :]
                    nearloci = nearloci[nearloci > -1]
                    prcp_max[r, c, :] = np.nanmax(prcpstn_ym[nearloci, :], axis=0)
                    # prcp_min[r, c, :] = np.nanmin(prcpstn_ym[nearloci, :], axis=0)
                # if near_grid_tempLoc[r, c, 0] > -1:
                #     nearloci = near_grid_tempLoc[r, c, :]
                #     nearloci = nearloci[nearloci > -1]
                #     tmean_max[r, c, :] = np.nanmax(tmeanstn_ym[nearloci, :], axis=0)
                #     trange_max[r, c, :] = np.nanmax(trangestn_ym[nearloci, :], axis=0)
                #     tmean_min[r, c, :] = np.nanmin(tmeanstn_ym[nearloci, :], axis=0)
                #     trange_min[r, c, :] = np.nanmin(trangestn_ym[nearloci, :], axis=0)
        prcp_max = au.transform(prcp_max, 3, 'box-cox')
        # prcp_min = au.transform(prcp_min, 3, 'box-cox')

        ################################################################################################################

        # calculate auto_corr and t_p_corr
        windows = 1  # parameters for auto-cc t-p-cc calculation: 1 could be better than 31
        lag = 1
        mean_autocorr, mean_tp_corr = \
            au.cc_calculate(windows, lag, prcpstn_ym, tmeanstn[:, indym], trangestn[:, indym])
        print('Tmean lag-1 daily autocorrelation: ', mean_autocorr)
        print('Trange-prcp daily correlation: ', mean_tp_corr)

        ################################################################################################################

        # # load OI-merged pop, pcp, tmean, trange
        fileoi = path_oi_pop + '/oimerge_pop' + str(year * 100 + month) + '.npz'
        datatemp = np.load(fileoi)
        pop = datatemp['oi_value']
        pop[pop<0] = 0
        pop[pop>1] = 1
        # pop_err = datatemp['oi_error']
        # temporal solution for pop
        # fileoi = '/datastore/GLOBALWATER/CommonData/EMDNA_new/ReanalysisCorrMerge/pop-old/bmamerge_pop_' + \
        #          str(year * 100 + month) + '.npz'
        # datatemp = np.load(fileoi)
        # pop = datatemp['bma_data']
        # # pop_err = datatemp['oi_error']

        fileoi = path_oi + '/oimerge_prcp' + str(year * 100 + month) + '_boxcox.npz'
        datatemp = np.load(fileoi)
        prcp = datatemp['oi_value']  # value in normal space
        prcp_err = datatemp['oi_error']

        fileoi = path_oi + '/oimerge_tmean' + str(year * 100 + month) + '.npz'
        datatemp = np.load(fileoi)
        tmean = datatemp['oi_value']
        tmean_err = datatemp['oi_error']

        fileoi = path_oi + '/oimerge_trange' + str(year * 100 + month) + '.npz'
        datatemp = np.load(fileoi)
        trange = datatemp['oi_value']
        trange[trange < 1] = 1
        trange_err = datatemp['oi_error']


        ################################################################################################################

        # constrain error value range because sometimes anomalous values will result in unrealistic errors
        # this is actually a very basic control with large value range allowed for error

        # lim1 = (prcp - prcp_max)
        # lim2 = (prcp - prcp_min)
        # ind = lim1 < lim2
        # lim = lim1
        # lim[ind] = lim2[ind]
        # ind = prcp_err > lim
        # prcp_err[ind] = lim[ind]
        #
        # lim1 = np.abs(tmean - tmean_max)
        # lim2 = np.abs(tmean - tmean_min)
        # ind = lim1 < lim2
        # lim = lim1
        # lim[ind] = lim2[ind]
        # ind = tmean_err > lim
        # tmean_err[ind] = lim[ind]
        #
        # lim1 = np.abs(trange - trange_max)
        # lim2 = np.abs(trange - trange_min)
        # ind = lim1 < lim2
        # lim = lim1
        # lim[ind] = lim2[ind]
        # ind = trange_err > lim
        # trange_err[ind] = lim[ind]

        ################################################################################################################

        # save to netcdf
        prcp_max = np.flipud(prcp_max)
        pop = np.flipud(pop)
        prcp = np.flipud(prcp)
        prcp_err = np.flipud(prcp_err)
        tmean = np.flipud(tmean)
        tmean_err = np.flipud(tmean_err)
        trange = np.flipud(trange)
        trange_err = np.flipud(trange_err)

        au.save_output_nc(fileout, gridinfo, seconds, mean_autocorr, mean_tp_corr, pop, prcp, tmean, trange,
                          prcp_err, tmean_err, trange_err, prcp_max)
