# prepare data for GMET probabilistic estimation based on OI merging and original station data
import numpy as np
import auxiliary as au
from auxiliary_merge import m_DateList
from calendar import monthrange
import datetime as dt
import netCDF4 as nc
import os
# control parameters
yearall = [1979,1979]
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
near_loc_grid = datatemp['near_grid_prcpLoc']
near_weight_grid = datatemp['near_grid_prcpWeight']
near_dist_grid = datatemp['near_grid_prcpDist']
near_loc_grid = np.flipud(near_loc_grid)
near_weight_grid = np.flipud(near_weight_grid)
near_dist_grid = np.flipud(near_dist_grid)

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

        # derive y_max for each grid pixels
        # load station data
        prcpstn_ym = prcpstn[:, indym]

        y_max = np.nan * np.zeros([nrows, ncols, ntimes], dtype=np.float32)
        for r in range(nrows):
            for c in range(ncols):
                if near_loc_grid[r, c, 0] < 0:
                    continue
                nearloci = near_loc_grid[r, c, :]
                nearloci = nearloci[nearloci > -1]
                y_max[r, c, :] = np.nanmax(prcpstn_ym[nearloci, :], axis=0)

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
        # fileoi = path_oi + '/oimerge_pop' + str(year * 100 + month) + '.npz'
        # datatemp = np.load(fileoi)
        # pop = datatemp['oi_value']
        # # pop_err = datatemp['oi_error']
        # temporal solution for pop
        fileoi = '/datastore/GLOBALWATER/CommonData/EMDNA_new/ReanalysisCorrMerge/pop-old/bmamerge_pop_' + \
                 str(year * 100 + month) + '.npz'
        datatemp = np.load(fileoi)
        pop = datatemp['bma_data']
        # pop_err = datatemp['oi_error']

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
        trange_err = datatemp['oi_error']

        ################################################################################################################

        # save to netcdf
        y_max = np.flipud(y_max)
        pop = np.flipud(pop)
        prcp = np.flipud(prcp)
        prcp_err = np.flipud(prcp_err)
        tmean = np.flipud(tmean)
        tmean_err = np.flipud(tmean_err)
        trange = np.flipud(trange)
        trange_err = np.flipud(trange_err)

        au.save_output_nc(fileout, gridinfo, seconds, mean_autocorr, mean_tp_corr, pop, prcp, tmean, trange,
                          prcp_err, tmean_err, trange_err, y_max)
