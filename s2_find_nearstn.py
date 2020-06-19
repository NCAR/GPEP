# for each target station and target grid cell, find its nearby stations
# calculation time: ~2 hours

import numpy as np
import netCDF4 as nc
import auxiliary as au
import os
import sys

# setting: searching nearby stations
nearstn_min = 20  # nearby stations: minimum number
nearstn_max = 30  # nearby stations: maximum number
search_radius = 400  # km. only search stations within this radius even nearstn_max cannot be reached
max_dist = 100  # max_distance in distance-based weight calculation

# input files
FileStnInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/stnlist_whole.txt'  # station basic information (lists)
FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/gridinfo_whole.nc'  # study area information
gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_original.npz'

# outfile
file_nearstn = '/home/gut428/nearstn_catalog.npz'

########################################################################################################################

print('Read station/grid information')
# read station information
stnID, stninfo = au.readstnlist(FileStnInfo)
nstn = len(stnID)

# read study area basic information
ncfid = nc.Dataset(FileGridInfo)
gridlat = ncfid.variables['latitude'][:].data
gridlon = ncfid.variables['longitude'][:].data
gridele = ncfid.variables['elev'][:].data
gridgns = ncfid.variables['gradient_n_s'][:].data
gridgwe = ncfid.variables['gradient_w_e'][:].data
mask = ncfid.variables['mask'][:].data  # 1: grids to be considered; the other values: invalid grids
ncfid.close()

nrows, ncols = np.shape(gridlat)
gridinfo = np.zeros([nrows, ncols, 6])
gridinfo[:, :, 0] = 1
gridinfo[:, :, 1] = gridlat
gridinfo[:, :, 2] = gridlon
gridinfo[:, :, 3] = gridele
gridinfo[:, :, 4] = gridgns
gridinfo[:, :, 5] = gridgwe
del gridlat, gridlon, gridele, gridgns, gridgwe

########################################################################################################################

# load station data
# station data are only used to indicate whether a station should be included in the catalog of nearby stations
if not os.path.isfile(gmet_stndatafile):
    print(gmet_stndatafile,'does not exist')
    sys.exit()
else:
    print('load station data')
    datatemp=np.load(gmet_stndatafile)
    prcp_stn_daily = datatemp['prcp_stn'][:, 0:2]
    tmean_stn_daily = datatemp['tmean_stn'][:, 0:2]
    trange_stn_daily = datatemp['trange_stn'][:, 0:2]

########################################################################################################################

# find nearby stations for stations/grids
near_grid_prcpLoc, near_grid_prcpDist, near_grid_prcpWeight, \
near_grid_tempLoc, near_grid_tempDist, near_grid_tempWeight, \
near_stn_prcpLoc, near_stn_prcpDist, near_stn_prcpWeight, \
near_stn_tempLoc, near_stn_tempDist, near_stn_tempWeight \
    = au.station_weight(prcp_stn_daily, tmean_stn_daily, stninfo, gridinfo, mask,
                        search_radius, nearstn_min, nearstn_max, max_dist)

# save data
np.savez_compressed(file_nearstn, near_grid_prcpLoc=near_grid_prcpLoc, near_grid_prcpDist=near_grid_prcpDist,
                    near_grid_prcpWeight=near_grid_prcpWeight, near_grid_tempLoc=near_grid_tempLoc,
                    near_grid_tempDist=near_grid_tempDist, near_grid_tempWeight=near_grid_tempWeight,
                    near_stn_prcpLoc=near_stn_prcpLoc, near_stn_prcpDist=near_stn_prcpDist,
                    near_stn_prcpWeight=near_stn_prcpWeight, near_stn_tempLoc=near_stn_tempLoc,
                    near_stn_tempDist=near_stn_tempDist, near_stn_tempWeight=near_stn_tempWeight)
