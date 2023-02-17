# the dynamic inputs in the test case are not compatible with xarray because of its lat/lon dims are the same with variables
# besides, it does not have a time axis. it forces every line corresponds to a time step, which is not flexible and will fail once start/end dates are changed
# this script makes a standard netcdf file based the test case

# will use ncrename

import subprocess
import xarray as xr
import numpy as np
import pandas as pd
import os, glob

inpath = './cali2017/griddata'
outpath = './cali2017/griddata_standard'
os.makedirs(outpath, exist_ok=True)

# change lat/lon
files = glob.glob(f'{inpath}/*.nc')
files.sort()
for f in files:
    fn = f.split('/')[-1]
    _ = subprocess.run(f'ncrename -v lat,lat2d -v lon,lon2d {f} {outpath}/{fn}', shell=True)


# merge
outfilelist = []
for i in range(4):
    infiles = []
    for j in range(i*7, (i+1)*7):
        infiles.append(f'{outpath}/201702{j+1:02}_HRRR_regrid_cube_pr.nc')
    dsi = xr.open_mfdataset(infiles, concat_dim='time', combine='nested')
    dsi['time'] = pd.date_range(f'2017-2-{i*7+1}', f'2017-2-{(i+1)*7}')
    dsi['lon'] = dsi['lon2d'].values[0][0,:]
    dsi['lat'] = dsi['lat2d'].values[0][:,0]
    dsi = dsi.drop_vars(['lon2d', 'lat2d'])
    outfile = f'{outpath}/inputfile_{i}.nc'
    dsi.to_netcdf(outfile)
    outfilelist.append(outfile)

with open(f'{outpath}/grid_file_list.txt', 'w') as f:
    for l in outfilelist:
        l = os.path.abspath(l)
        _ = f.write(l+'\n')
