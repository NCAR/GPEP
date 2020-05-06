import numpy as np
import netCDF4 as nc
import auxiliary as au
import regression as reg
import datetime as dt
from matplotlib import pyplot as plt
from scipy import io
import os
import sys
from scipy.interpolate import griddata


def ncread(file, var):
    # read a variable from netcdf
    ncfid = nc.Dataset(file)
    data = ncfid[var][:].data
    ncfid.close()
    return data


def neargrid(rowtar, coltar, rowori, colori, hwsize):
    # inputs are 1D matrices
    # tar is target area
    # ori is original area
    # hwsize is half window size (e.g., 4 means the space window width/length is 2*4+1)
    # find a space window centering the target grid in the original area and calculate the weights
    nrows = len(rowtar)
    ncols = len(coltar)
    rowse = np.nan * np.zeros([nrows, ncols, 2])  # se: start/end
    colse = np.nan * np.zeros([nrows, ncols, 2])  # se: start/end
    weight = np.nan * np.zeros([nrows, ncols, (hwsize * 2 + 1) ** 2])  # from left to right/from top to bottom weight

    for rr in range(nrows):
        rowloc = np.argmin(np.abs(rowori - rowtar[rr]))
        rowse[rr, :, 0] = rowloc - hwsize
        rowse[rr, :, 1] = rowloc + hwsize

    for cc in range(ncols):
        colloc = np.argmin(np.abs(colori - coltar[cc]))
        colse[:, cc, 0] = colloc - hwsize
        colse[:, cc, 1] = colloc + hwsize

    rowse[rowse < 0] = 0
    rowse[rowse > nrows] = nrows
    colse[colse < 0] = 0
    colse[colse > ncols] = nrows

    maxdist = (hwsize + 0.5) * np.sqrt(2) + 0.5
    for rr in range(nrows):
        for cc in range(ncols):
            rowse_rc = rowse[rr, cc, :]
            colse_rc = colse[rr, cc, :]
            flag = 0
            for i in range(rowse_rc[0], rowse_rc[1] + 1):
                for j in range(colse_rc[0], colse_rc[1] + 1):
                    dist = ((i - rr) ** 2 + (j - cc) ** 2) ** 0.5
                    weight[rr, cc, flag] = au.distanceweight(dist, maxdist, 3)
                    flag = flag + 1

    return rowse, colse, weight


# basic information
year = [2018, 2018]
month = [1,12]
lontar = np.arange(-180 + 0.05, -50, 0.1)
lattar = np.arange(5 + 0.05, 80, 0.1)
hwsize = 2

# ERA-5
inpath = ''
outpath = ''
filenear = outpath + '/neargrid.npz'
var = 'tp'
for y in range(year[0], year[1] + 1):
    for m in range(month[0],month[1]+1):
        file = inpath + '/ERA5_' + str(y) + '.nc'
        # read location information
        if not os.path.isfile(filenear):
            latori = ncread(file,'latitude')
            lonori = ncread(file, 'longitude')
            rowse, colse, weight = neargrid(lattar,lontar,latori,lonori,hwsize)
            io.savemat(filenear,{'rowse':rowse, 'colse':colse, 'weight':weight})
        else:
            datatemp = io.loadmat(filenear)
            rowse = datatemp['rowse']
            colse = datatemp['colse']
            weight = datatemp['weight']
        # read original ERA5 data and sum it to daily scale
        
