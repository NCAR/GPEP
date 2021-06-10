import xarray as xr
import pandas as pd
import glob, sys, os, time

def read_EMDNA(file):
    with xr.open_dataset(file) as ds0:
        ds = ds0.copy()
        ds.coords['x'] = ds.longitude
        ds.coords['y'] = ds.latitude
        year = int(ds.date[0]/10000)
        ds.coords['time'] = pd.date_range(start=f'{year}-01-01', end=f'{year}-12-31', freq='1D')
    return ds

year = int(sys.argv[1])
# year = 2000
path = f'/project/rpp-kshook/gut428/EMDNA_v1/Estimate/{year}'

filelist = glob.glob(path+'/*.nc4')
for file_in in filelist:
    file_out = file_in.replace('nc4', 'nc')
    if os.path.isfile(file_out):
        continue
    t1 = time.time()
    print('----------------------------------------------------------------')
    print('infile:', file_in)
    print('outfile:', file_out)
    ds = read_EMDNA(file_in)
    ds = ds.drop_vars(['latitude','longitude'])
    encoding = {}
    for key in ds.keys():
        encoding[key] = {'zlib': True, 'complevel': 5, '_FillValue': -9999.0}
    ds.to_netcdf(file_out, encoding=encoding)
    t2 = time.time()
    print(t2-t1)