# prepare data for GMET probabilistic estimation based on OI merging and original station data
import numpy as np
import auxiliary as au
from auxiliary_merge import m_DateList
from calendar import monthrange
import datetime as dt
import netCDF4 as nc

# control parameters
year = 2018
month = 2

### Mac settings
path_oi = '/Users/localuser/Research/EMDNA/oimerge'
near_file_GMET = '/Users/localuser/Research/EMDNA/regression/weight_nearstn.npz'
gmet_stndatafile = '/Users/localuser/Research/EMDNA/stndata_whole.npz'
FileStnInfo = '/Users/localuser/Research/EMDNA/basicinfo/stnlist_whole.txt'  # station basic information (lists)
FileGridInfo = '/Users/localuser/Research/EMDNA/basicinfo/gridinfo_whole.nc'  # study area information
### Mac settings

# ### Plato settings
# path_oi = '/home/gut428/OImerge'
# near_file_GMET = '/datastore/GLOBALWATER/CommonData/EMDNA/PyGMETout/weight.npz'
# gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA/stndata_whole.npz'
# FileStnInfo = '/datastore/GLOBALWATER/CommonData/EMDNA/StnGridInfo/stnlist_whole.txt'  # station basic information (lists)
# FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA/StnGridInfo/gridinfo_whole.nc'  # study area information
# ### Plato settings

# output files
FileRegression = '/home/gut428/GMET_regression' + '/reg_' + str(year*100+month) + '.nc'

########################################################################################################################

# 1. basic information

print('Read study area basic information')
# station location and attribute information
# stninfo: [ stations, 1/lat/lon/elev/slope_ns/slope_we ]
stnID, stninfo = au.readstnlist(FileStnInfo)
nstn = len(stnID)

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

# grid information
lontar = np.arange(-180 + 0.05, -50, 0.1)
lattar = np.arange(85 - 0.05, 5, -0.1)
nrows = len(lattar)
ncols = len(lontar)

########################################################################################################################

# 2. read study area basic information
print('Read study area basic information')
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

########################################################################################################################

# load near station information
datatemp = np.load(near_file_GMET)
near_loc_stn = datatemp['near_stn_prcpLoc']
near_weight_stn = datatemp['near_stn_prcpWeight']
near_dist_stn = datatemp['near_stn_prcpDist']
near_loc_grid = datatemp['near_grid_prcpLoc']
near_weight_grid = datatemp['near_grid_prcpWeight']
near_dist_grid = datatemp['near_grid_prcpDist']
near_loc_grid = np.flipud(near_loc_grid)
near_weight_grid = np.flipud(near_weight_grid)
near_dist_grid = np.flipud(near_dist_grid)

########################################################################################################################

# derive y_max for each grid pixels

# load station data
date_list, date_number = m_DateList(1979, 2018, 'ByYear')
indym = (date_number['yyyy'] == year) & (date_number['mm'] == month)
datatemp = np.load(gmet_stndatafile)
stndata = datatemp['prcp_stn'][:, indym]
stnlle = datatemp['stn_lle']
nstn, ntimes = np.shape(stndata)
del datatemp

y_max = np.zeros([nrows, ncols, ntimes], dtype=np.float32)
for r in range(nrows):
    for c in range(ncols):
        if near_loc_grid[r, c, 0] < 0:
            continue
        nearloci = near_loc_grid[r, c, :]
        nearloci = nearloci[nearloci > -1]
        y_max[r, c, :] = np.nanmax(stndata[nearloci, :], axis=0)

########################################################################################################################

# calculate auto_corr and t_p_corr
windows = 1  # parameters for auto-cc t-p-cc calculation: 1 could be better than 31
lag = 1
datatemp = np.load(gmet_stndatafile)
stndata = datatemp['prcp_stn'][:, indym]
mean_autocorr, mean_tp_corr = \
    au.cc_calculate(windows, lag, datatemp['prcp_stn'][:, indym],
                    datatemp['tmean_stn'][:, indym], datatemp['trange_stn'][:, indym])
print('Tmean lag-1 daily autocorrelation: ', mean_autocorr)
print('Trange-prcp daily correlation: ', mean_tp_corr)

########################################################################################################################
# load OI-merged pop, pcp, tmean, trange
fileoi = path_oi + '/oimerge_pop' + str(year * 100 + month) + '.npz'
datatemp = np.load(fileoi)
pop = datatemp['oi_value']
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


# flip data to accommodate grid information
pop = np.flipud(pop)
prcp = np.flipud(prcp)
prcp_err = np.flipud(prcp_err)
tmean = np.flipud(tmean)
tmean_err = np.flipud(tmean_err)
trange = np.flipud(trange)
trange_err = np.flipud(trange_err)

########################################################################################################################
# save to netcdf
au.save_output_nc(FileRegression, gridinfo, seconds, mean_autocorr, mean_tp_corr, pop, prcp, tmean, trange,
                  prcp_err, tmean_err, trange_err, y_max)
