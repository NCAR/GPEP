# Get a test case from Fortran GMET repo
# https://github.com/NCAR/GMET/tree/master/test_cases

import subprocess
import xarray as xr
import pandas as pd
import os, glob

def prepare_python_input():

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
        for j in range(i * 7, (i + 1) * 7):
            infiles.append(f'{outpath}/201702{j + 1:02}_HRRR_regrid_cube_pr.nc')
        dsi = xr.open_mfdataset(infiles, concat_dim='time', combine='nested')
        dsi['time'] = pd.date_range(f'2017-2-{i * 7 + 1}', f'2017-2-{(i + 1) * 7}')
        dsi['lon'] = dsi['lon2d'].values[0][0, :]
        dsi['lat'] = dsi['lat2d'].values[0][:, 0]
        dsi = dsi.drop_vars(['lon2d', 'lat2d'])
        outfile = f'{outpath}/inputfile_{i}.nc'
        dsi.to_netcdf(outfile)
        outfilelist.append(outfile)

    with open(f'{outpath}/grid_file_list.txt', 'w') as f:
        for l in outfilelist:
            l = os.path.abspath(l)
            _ = f.write(l + '\n')


targetpath = '.'
os.chdir(targetpath)

# get test case and delete unnecessary files
os.system('git clone git@github.com:NCAR/GMET.git')
os.system('tar -xf  GMET/test_cases/cali2017.tgz')
os.system('rm -rf GMET')

# process fortran gridded inputs (e.g., add time axis)
prepare_python_input()

