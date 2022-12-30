# all data processing related functions

import os, time
import pandas as pd
import numpy as np
import xarray as xr

def assemble_fortran_GMET_stns_to_one_file(config):
    # Fortran GMET assumes that each station has an independent file. PyGMET assumes that point-based stations are
    # small concerning data size, and thus it is more convenient to assemble all stations in one file
    t1 = time.time()

    # parse and change configurations
    outpath_parent = config['outpath_parent']
    path_stn_info = f'{outpath_parent}/stn_info'
    os.makedirs(path_stn_info, exist_ok=True)
    file_allstn = f'{path_stn_info}/all_stn.nc'

    config['path_stn_info'] = path_stn_info
    config['file_allstn'] = file_allstn

    # in/out information to this function
    input_stn_list = config['input_stn_list']
    input_stn_path = config['input_stn_path']
    file_allstn = config['file_allstn']

    # settings and prints
    minTrange = 1  # hard-coded setting
    var_names = ['prcp', 'tmin', 'tmax'] # variables in the files
    add_tmean_trange = True # if tmin and tmax are present, calculate tmean and trange

    print('#' * 50)
    print('Assembling individual station files to one single file')
    print('#' * 50)
    print('Input station list:', input_stn_list)
    print('Input station folder:', input_stn_path)
    print('Output station file:', file_allstn)
    print('Target variables:', var_names)
    print('add_tmean_trange:', add_tmean_trange)
    if add_tmean_trange:
        print('minTrange:', minTrange)

    if os.path.isfile(file_allstn):
        print('Note! Output station file exists')
        if config['overwrite_stninfo'] == True:
            print('overwrite_stninfo is True. Continue.')
        else:
            print('overwrite_stninfo is False. Skip station assembling.')
            return config

    # start assembling
    df_stn = pd.read_csv(input_stn_list)
    nstn = len(df_stn)

    infile0 = f'{input_stn_path}/{df_stn.stnid.values[0]}.nc'
    with xr.open_dataset(infile0) as ds:
        time_coord = ds.time
        ntime = len(ds.time)

    var_values = np.nan * np.zeros([nstn, ntime, len(var_names)], dtype=np.float32)
    print(f'Station number: {nstn}')
    print(f'Time steps: {ntime}')
    print('Start assembling ... ')

    for i in range(nstn):
        infilei = f'{input_stn_path}/{df_stn.stnid.values[i]}.nc'
        dsi = xr.load_dataset(infilei)
        for j in range(len(var_names)):
            if var_names[j] in dsi.data_vars:
                var_values[i, :, j] = np.squeeze(dsi[var_names[j]].values)

    ds_stn = xr.Dataset()
    ds_stn.coords['time'] = time_coord
    ds_stn.coords['stn'] = np.arange(nstn)

    for col in df_stn.columns:
        ds_stn[col] = xr.DataArray(df_stn[col].values, dims=('stn'))

    for i in range(len(var_names)):
        ds_stn[var_names[i]] = xr.DataArray(var_values[:, :, i], dims=('stn', 'time'))

    # add tmean and trange
    if (add_tmean_trange == True) and ('tmin' in ds_stn.data_vars) and ('tmax' in ds_stn.data_vars):
        ds_stn['tmean'] = (ds_stn['tmin'] + ds_stn['tmax']) / 2
        ds_stn['trange'] = np.abs(ds_stn['tmax'] - ds_stn['tmin'])
        # constrain trange
        v = ds_stn['trange'].values
        v[v < minTrange] = minTrange
        ds_stn['trange'].values = v

    # save to output files
    encoding = {}
    for var in ds_stn.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}
    ds_stn.to_netcdf(file_allstn, encoding=encoding)

    t2 = time.time()
    print('Time cost (seconds):', t2-t1)
    print('Successful assembling!')

    return config


