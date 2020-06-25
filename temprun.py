# optimal interpolation merging
# merge background (reanalysis) and observation (regression estimates)
# computation time:
# oi-merge for stations: 3 variables adn 12 months = 36 jobs. ~0.5 hour per job
# oi-merge for grids: 3 variables adn 12 months = 36 jobs. ~0.5 hour per job

import numpy as np
from scipy import io
import os
import netCDF4 as nc
from optimal_interpolation import OImerge
from auxiliary_merge import *

########################################################################################################################

# time periods and methods
vars = sys.argv[1]
vars = [vars]
month = int(sys.argv[2])

print(vars,month)

########################################################################################################################

# basic settings
weightmode = 'BMA' # method used to merge different reanalysis products
# vars = ['prcp', 'tmean', 'trange']

### Plato settings
# input files/paths
FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/gridinfo_whole.nc'
FileStnInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/stnlist_whole.txt'
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_aftercheck.npz'
path_bac = '/datastore/GLOBALWATER/CommonData/EMDNA_new/ReanalysisCorrMerge/pop'
path_obs = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck'
near_file_GMET = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck/nearstn_catalog.npz'
file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA_new/DEM/NA_DEM_010deg_trim.mat'

# output files/paths (can also be used as inputs once generated)
path_oimerge = '/home/gut428/OImerge_GWRLSBMA'
### Plato settings

file_regression_stn = path_obs + '/regression_stn.npz'
file_corrmerge_stn = [''] * len(vars)
for i in range(len(vars)):
    if vars[i] == 'pop':
        file_corrmerge_stn[i] = path_bac + '/merge_stn_' + vars[i] + '_GWR_' + weightmode + '.npz'
    else:
        file_corrmerge_stn[i] = path_bac + '/mergecorr_stn_' + vars[i] + '_GWRLS_' + weightmode + '.npz'

########################################################################################################################

# basic processing
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

########################################################################################################################

# OI-merging at grid scale

for v in range(len(vars)):
    print('OI merge at grids:', vars[v])

    # load station original observations
    datatemp = np.load(gmet_stndatafile)
    if vars[v] == 'pop':
        observation_stn = datatemp['prcp_stn']
        observation_stn[observation_stn > 0] = 1
    else:
        observation_stn = datatemp[vars[v] + '_stn']
    del datatemp

    # load near station information
    datatemp = np.load(near_file_GMET)
    if vars[v] == 'prcp' or vars[v] == 'pop':
        near_loc = datatemp['near_grid_prcpLoc']
        near_weight = datatemp['near_grid_prcpWeight']
        near_dist = datatemp['near_grid_prcpDist']
    else:
        near_loc = datatemp['near_grid_tempLoc']
        near_weight = datatemp['near_grid_tempWeight']
        near_dist = datatemp['near_grid_tempDist']
    near_loc = np.flipud(near_loc)
    near_weight = np.flipud(near_weight)
    near_dist = np.flipud(near_dist)
    del datatemp

    # start OI merging
    for m in range(month-1, month):
        print('month', m + 1)
        indm = (date_mm == m + 1)
        nday = sum(indm)
        datem = date_yyyy[indm]

        # load OI merged data at station points
        filemerge_stn = '/datastore/GLOBALWATER/CommonData/EMDNA_new/OImerge_GWRLSBMA/OImerge_stn_GWRLSBMA_' + vars[v] + '.npz'
        datatemp = np.load(filemerge_stn)
        oimerge_stn = datatemp['oimerge_stn']
        del datatemp

        for y in range(1979, 2019):
            print('year',y)
            fileoi_ym_exist = '/datastore/GLOBALWATER/CommonData/EMDNA_new/OImerge_GWRLSBMA/oimerge_' + vars[v] + str(y * 100 + m + 1) + '.npz'
            fileoi_ym = path_oimerge + '/oimerge_' + vars[v] + str(y*100+m+1) + '.npz'
            indym1 = datem == y
            ndayy = np.sum(indym1)
            indym2 = (date_mm == m + 1) & (date_yyyy == y)

            if vars[v] == 'prcp':
                # value and error in normal space
                fileoi_ym_boxcox = path_oimerge + '/oimerge_' + vars[v] + str(y * 100 + m + 1) + '_boxcox.npz'
                if os.path.isfile(fileoi_ym_boxcox):
                    continue
                transmode = 'box-cox'
                tranexp = 4
                oi_error_stn = ( au.transform(oimerge_stn[:, indym2],tranexp,transmode) -
                                au.transform(observation_stn[:, indym2], tranexp, transmode) ) ** 2
                oi_error_grid = extrapolation(oi_error_stn, near_loc, near_dist, excflag=1)
                oi_error_grid = oi_error_grid ** 0.5
                oi_value = au.transform(oi_value, tranexp, transmode)
                np.savez_compressed(fileoi_ym_boxcox, oi_value=oi_value, oi_error=oi_error_grid,
                                    latitude=lattar, longitude=lontar, tranexp=tranexp, transmode=transmode)