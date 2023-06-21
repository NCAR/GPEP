# Get a test case from Fortran GMET repo (2023-3-21), do preprocessing, and delete unnecessary files
#   optional step -- independent test cases are included in the test_cases folder

import subprocess
import xarray as xr
import pandas as pd
import os, glob, pathlib, shutil

def process_GMET_dynamic_inputs(inpath, outpath):

    os.makedirs(outpath, exist_ok=True)

    # change lat/lon
    files = glob.glob(f'{inpath}/*.nc')
    files.sort()
    for f in files:
        fn = f.split('/')[-1]
        _ = subprocess.run(f'ncrename -v lat,lat2d -v lon,lon2d {f} {outpath}/{fn}', shell=True)

    # merge
    outfilelist = []
    filenum = 2
    timenum = 5
    for i in range(filenum):
        infiles = []
        for j in range(i * timenum, (i + 1) * timenum):
            infiles.append(f'{outpath}/201702{j + 1:02}_HRRR_regrid_cube_pr.nc')
        dsi = xr.open_mfdataset(infiles, concat_dim='time', combine='nested')
        dsi['time'] = pd.date_range(f'2017-2-{i * timenum + 1}', f'2017-2-{(i + 1) * timenum}')
        dsi['lon'] = dsi['lon2d'].values[0][0, :]
        dsi['lat'] = dsi['lat2d'].values[0][:, 0]
        dsi = dsi.drop_vars(['lon2d', 'lat2d'])
        outfile = f'{outpath}/inputfile_{i}.nc'
        dsi.to_netcdf(outfile)
        outfilelist.append(outfile)

    _ = os.system(f'rm {outpath}/201702*_HRRR_regrid_cube_pr.nc')

    # generate a list file
    outpath0 = str(pathlib.Path(outpath).parent)
    with open(f'{outpath0}/gridded_datafile_list.txt', 'w') as f:
        for l in outfilelist:
            # l = os.path.abspath(l)
            l = pathlib.Path(l).name
            _ = f.write(f'./griddata/{l}\n')


def reduce_stn_data_length(inpath):
    files = glob.glob(f'{inpath}/*.nc')
    files.sort()
    for f in files:
        ds = xr.load_dataset(f)
        ds = ds.isel(time=slice(0, 10))
        ds.to_netcdf(f+'2')
        _ = os.remove(f)
        _ = os.system(f'mv {f}2 {f}')

def prepare_python_input():

    # gridded data (i.e., dynamic predictors)
    inpath = './cali2017/griddata'
    outpath = './cali2017/inputs/griddata'
    process_GMET_dynamic_inputs(inpath, outpath)

    # # station data
    # reduce_stn_data_length('./cali2017/stndata')

    # remove redundant files
    _ = shutil.rmtree('./cali2017/griddata')
    _ = shutil.rmtree('./cali2017/outputs')
    _ = shutil.rmtree('./cali2017/run')
    _ = os.remove('./cali2017/inputs/CALI.screened_stn_list_slope.v3.csv.prev.gz')
    _ = os.remove('./cali2017/inputs/hrrr_cube_prec_list.txt')


    # move folders
    _ = os.system(f'mv ./cali2017/stndata ./cali2017/inputs')

targetpath = '../test_cases/'
os.chdir(targetpath)

# get test case and delete unnecessary files
os.system('git clone git@github.com:NCAR/GMET.git')
os.system('tar -xf  GMET/test_cases/cali2017.tgz')
os.system('rm -rf GMET')

# process fortran gridded inputs (e.g., add time axis)
prepare_python_input()

