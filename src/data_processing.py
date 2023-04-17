# all data processing related functions

import os, time, sys
import pandas as pd
import numpy as np
import xarray as xr

########################################################################################################################
# data transformation

def boxcox_transform(data, texp=4):
    # transform prcp to approximate normal distribution
    # mode: box-cox; power-law
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    data[data<0] = 0
    datat = (data ** (1 / texp) - 1) / (1 / texp)
    # datat[data < -3] = -3
    return datat

def boxcox_back_transform(data, texp=4):
    # transform prcp to approximate normal distribution
    # mode: box-cox; power-law
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    data[data<-texp] = -texp
    datat = (data / texp + 1) ** texp
    return datat

def data_transformation(data, method, settings, mode='transform'):
    if method == 'boxcox':
        if mode == 'transform':
            data = boxcox_transform(data, settings['exponent'])
        elif mode == 'back_transform':
            data = boxcox_back_transform(data, settings['exponent'])
        else:
            print('Unknown transformation mode: entry=', mode); sys.exit()
    else:
        print('Unknown transformation method: entry=', method); sys.exit()
    return data


########################################################################################################################
# input station data processing

def merge_stndata_into_single_file(config):
    # GMET v2.0 assumes that each station has an independent file. PyGMET will merge the station data into one file,
    #   which can speed up i/o in subsequent runs of PYGMET using the same dataset. 
    t1 = time.time()

    # parse and change configurations
    outpath_parent = config['outpath_parent']
    path_stn_info  = f'{outpath_parent}/stn_info'
    file_allstn    = f'{path_stn_info}/all_station_data.nc'     # set default name for a merged station data file
    os.makedirs(path_stn_info, exist_ok=True)

    config['path_stn_info'] = path_stn_info
    config['file_allstn']   = file_allstn             # store default name for a merged station data file in config

    # in/out information to this function
    if 'input_stn_list' in config:
        input_stn_list = config['input_stn_list']
    else:
        input_stn_list = ''
    if 'input_stn_path' in config:
        input_stn_path = config['input_stn_path']
    else:
        input_stn_path = ''
    if 'input_stn_all' in config:
        input_stn_all  = config['input_stn_all']
    else:
        input_stn_all  = ''

    file_allstn = config['file_allstn']
    input_vars  = config['input_vars']
    target_vars = config['target_vars']

    if 'minRange_vars' in config:
        minRange_vars = config['minRange_vars']
        if isinstance(minRange_vars, (int, float)):
            minRange_vars = [minRange_vars] * len(target_vars)
    else:
        minRange_vars = [-np.inf] * len(target_vars)

    if 'maxRange_vars' in config:
        maxRange_vars = config['maxRange_vars']
        if isinstance(maxRange_vars, (int, float)):
            maxRange_vars = [maxRange_vars] * len(target_vars)
    else:
        minRange_vars = [np.inf] * len(target_vars)

    if 'transform_vars' in config:
        transform_vars = config['transform_vars']
        if not isinstance(transform_vars, list):
            transform_vars = [transform_vars] * len(target_vars)
    else:
        transform_vars = [''] * len(target_vars)

    transform_settings = config['transform']
    mapping_InOut_var  = config['mapping_InOut_var']

    if 'overwrite_stninfo' in config:
        overwrite_stninfo = config['overwrite_stninfo']
    else:
        overwrite_stninfo = False
        
    if 'overwrite_merged_stnfile' in config:
        overwrite_merged_stnfile = config['overwrite_merged_stnfile']
    else:
        overwrite_merged_stnfile = True

    # settings and prints
    print('#' * 50)
    print('Merging individual station files to one single file')
    print('#' * 50)
    print('Input station list:  ', input_stn_list)
    print('Input station folder:', input_stn_path)
    print('Output station file: ', file_allstn)
    print('Target variables:    ', input_vars)

    if os.path.isfile(file_allstn):
        print('NOTE: Merged station file exists')
        if overwrite_merged_stnfile == True:
            print('overwrite_merged_stnfile is True. Continue.')
        else:
            print('overwrite_merged_stnfile is False. Skip station merging')
            return config

    # load station data
    if 'input_stn_all' in config and os.path.isfile(input_stn_all):
        print('input_stn_all exists:    ', input_stn_all)
        print('reading station info from', input_stn_all, 'instead of individual files')
        #print('reading station information from', file_allstn, 'instead of individual files')  # this looks incorrect

        ds_stn = xr.load_dataset(input_stn_all)

        if os.path.isfile(input_stn_list):
            df_stn = pd.read_csv(input_stn_list)
            for col in df_stn.columns:
                ds_stn[col] = xr.DataArray(df_stn[col].values, dims=('stn'))

    else:
        print('A merged station data file does not exist: reading station information from individual files')

        df_stn = pd.read_csv(input_stn_list)
        nstn   = len(df_stn)

        infile0 = f'{input_stn_path}/{df_stn.stnid.values[0]}.nc'
        with xr.open_dataset(infile0) as ds:
            time_coord = ds.time
            ntime = len(ds.time)

        var_values = np.nan * np.zeros([nstn, ntime, len(input_vars)], dtype=np.float32)
        print(f'Number of input station: {nstn}')
        print(f'Time steps:              {ntime}')
        print('Start merging station data ... ')

        for i in range(nstn):
            infilei = f'{input_stn_path}/{df_stn.stnid.values[i]}.nc'
            dsi = xr.load_dataset(infilei)
            for j in range(len(input_vars)):
                if input_vars[j] in dsi.data_vars:
                    var_values[i, :, j] = np.squeeze(dsi[input_vars[j]].values)

        ds_stn                = xr.Dataset()
        ds_stn.coords['time'] = time_coord
        ds_stn.coords['stn']  = np.arange(nstn)

        for col in df_stn.columns:
            ds_stn[col] = xr.DataArray(df_stn[col].values, dims=('stn'))

        for i in range(len(input_vars)):
            ds_stn[input_vars[i]] = xr.DataArray(var_values[:, :, i], dims=('stn', 'time'))

        del var_values

    # convert input vars to target vars
    for mapping in mapping_InOut_var:
        mapping = mapping.replace(' ', '')
        tarvar  = mapping.split('=')[0]

        if tarvar in target_vars:

            if tarvar in ds_stn:
                print(f'{tarvar} is already in ds_stn. no need to perform {mapping}')
            else:
                invars = [v for v in input_vars if v in mapping]
                for v in invars:
                    exec(f"{v}=ds_stn.{v}")

                operation      = mapping.split('=')[1]
                ds_stn[tarvar] = eval(operation)

    # constrain variables
    for i in range(len(target_vars)):
        vari = target_vars[i]
        if vari in ds_stn.data_vars:
            v = ds_stn[vari].values
            if np.any(v < minRange_vars[i]):
                print(f'{vari} station data has values < {minRange_vars[i]}. Adjust those to {minRange_vars[i]}.')
                v[v < minRange_vars[i]] = minRange_vars[i]
            if np.any(v > maxRange_vars[i]):
                print(f'{vari} station data has values > {maxRange_vars[i]}. Adjust those to {maxRange_vars[i]}.')
                v[v > maxRange_vars[i]] = maxRange_vars[i]
            ds_stn[vari].values = v

    # transform variables
    print('Transform variables if relevant settings are provided')
    for i in range(len(transform_vars)):
        if len(transform_vars[i])>0:
            tvar = target_vars[i] + '_' + transform_vars[i]
            print(f'Perform {transform_vars[i]} transformation for {target_vars[i]}. Add a new variable {tvar} to output station file.')
            if tvar in ds_stn:
                print(f'{tvar} exists in ds_stn. no need to perform transformation')
                continue
            ds_stn[tvar] = ds_stn[vari].copy()
            ds_stn[tvar].values = data_transformation(ds_stn[target_vars[i]].values, transform_vars[i],
                                                      transform_settings[transform_vars[i]], 'transform')
        else:
            print(f'Do not perform transformation for {target_vars[i]}')

    # save to output files
    encoding = {}
    for var in ds_stn.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}

    ds_stn['stnid'] = ds_stn['stnid'].astype('|S')
    ds_stn.to_netcdf(file_allstn, encoding=encoding)

    t2 = time.time()
    print('Time cost (s) for merging station data file:', t2-t1, '\n')

    return config


