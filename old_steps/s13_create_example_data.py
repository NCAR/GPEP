# the complete EMDNA dataset is too large
# this script creates example files for a small region

import xarray as xr
import numpy as np

path='/project/6008034/Model_Output/ClimateForcingData/EMDNA.subset.CH/2015/'
outfile = 'EMDNA_example_2015.nc'
year = 2015

ens = [1, 25]
for e in range(ens[0], ens[1]+1):
    filee = f'{path}/EMDNA_{year}.{e:03}.CH.nc4'
    with xr.open_dataset(filee) as ds:
        ds = ds.expand_dims(dim='ens')
        if e == ens[0]:
            ds_out = ds
        else:
            ds_out = xr.concat((ds_out, ds), dim='ens')

ds_out.coords['ens'] = np.arange(ens[0], ens[1]+1)
ds_out = ds_out.rename({'x':'lon', 'y':'lat'})
comp = dict(zlib=True, complevel=5)
encoding = {var: comp for var in ds_out.data_vars}
ds_out.to_netcdf(outfile, encoding=encoding)
