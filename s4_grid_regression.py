# grid regression like GMET

import numpy as np
import regression as reg
from scipy import io
import os, sys
from bma_merge import bma
from auxiliary_merge import extrapolation
import netCDF4 as nc

FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/gridinfo_whole.nc'
gmet_stnfile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_aftercheck.npz'
file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA_new/DEM/NA_DEM_010deg_trim.mat'
near_file_GMET = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck/nearstn_catalog.npz'


########################################################################################################################

# basic processing
print('start basic processing')

# mask
mask = io.loadmat(file_mask)
mask = mask['DEM']
mask[~np.isnan(mask)] = 1  # 1: valid pixels
nrows, ncols = np.shape(mask)

# meshed lat/lon of the target region
ncfid = nc.Dataset(FileGridInfo)
lattarm = ncfid.variables['latitude'][:].data
lattarm = np.flipud(lattarm)
lontarm = ncfid.variables['longitude'][:].data
ncfid.close()
lontarm[np.isnan(mask)] = np.nan
lattarm[np.isnan(mask)] = np.nan
lontar = lontarm[0, :]
lattar = lattarm[:, 0]

# load observations for all stations
datatemp = np.load(gmet_stndatafile)
stndata = datatemp['prcp_stn']
stninfo = datatemp['stninfo']
stnID = datatemp['stnID']
date_ymd = datatemp['date_ymd']
nstn, ntimes = np.shape(stndata)
del datatemp
date_yyyy = (date_ymd/10000).astype(int)
date_mm = (np.mod(date_ymd, 10000)/100).astype(int)

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
del datatemp



