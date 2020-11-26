# subset the EMDNA according to lat/lon extent
import xarray as xr
import numpy as np
import pandas as pd
import os, sys

def generate_filelist(path, year, ens, suffix):
    num = ens[1] - ens[0] + 3 # include _mean and _spread
    filelist = [' '] * num
    for e in range(ens[0], ens[1]+1):
        filelist[e-1] = f'{path}/EMDNA_{year}.{e:03}{suffix}.nc4'
    filelist[-2] = f'{path}/EMDNA_{year}_mean{suffix}.nc4'
    filelist[-1] = f'{path}/EMDNA_{year}_spread{suffix}.nc4'
    return filelist

def read_EMDNA(file):
    ds = xr.open_dataset(file)
    ds.coords['x'] = ds.longitude
    ds.coords['y'] = ds.latitude
    year = int(ds.date[0]/10000)
    ds.coords['time'] = pd.date_range(start=f'{year}-01-01', end=f'{year}-12-31', freq='1D')
    ds.close()
    return ds


# parent path of EMDNA
emdna_path = '/home/gut428/projects/rpp-kshook/gut428/EMDNA_v1/Estimate'

# outpath for saving the subset data
suffix = '.CH'
subset_path = '/home/gut428/scratch/' + 'subset' + suffix

# parameters
yearrange = [1979, 2018] # should be within 【1979, 2019]
latrange = [50.7, 51.9] # should be within [5, 85]
lonrange = [-116.8, -115.2] # should be within [-180, -50]
ensrange = [1, 100] # should be within [1, 100]

# start
if not os.path.isdir(subset_path):
    os.mkdir(subset_path)

# just for saving output
encoding = {'prcp': {"zlib": True, "complevel": 5},
            'tmean': {"zlib": True, "complevel": 5},
            'trange': {"zlib": True, "complevel": 5}}

for year in range(yearrange[0], yearrange[1]+1):
    # input file list
    pathy_in = emdna_path + '/' + str(year)
    filelist_in = generate_filelist(pathy_in, year, ensrange, '')
    # output file list
    pathy_out = subset_path + '/' + str(year)
    if not os.path.isdir(pathy_out):
        os.mkdir(pathy_out)
    filelist_out = generate_filelist(pathy_out, year, ensrange, suffix)
    for filein, fileout in zip(filelist_in, filelist_out):
        # print information
        print('-'*100)
        print('InFile:',filein)
        print('OutFile:', fileout)
        if os.path.isfile(fileout):
            print('Outfile exists. Continue ...')
            continue
        # read raw EMDNA data
        datain = read_EMDNA(filein)
        # extract subset region
        dataout = datain.sel(y=slice(latrange[1], latrange[0]), x=slice(lonrange[0], lonrange[1]))
        # write data
        dataout = dataout[['prcp','tmean','trange']]
        dataout.to_netcdf(fileout, encoding=encoding)
        del dataout, datain




