# all data processing related functions

import os, time, sys
import pandas as pd
import numpy as np
import xarray as xr
from decorators import timer
from scipy.stats import norm
from scipy.interpolate import interp1d
from empirical_cdf import calculate_monthly_cdfs, normal_quantile_transform, inverse_normal_quantile_transform

########################################################################################################################
# data transformation

def boxcox_transform(data, texp=4):
    """
    Apply the Box-Cox transformation to the input data.
    
    Parameters:
    - data (array-like): Input data to be transformed.
    - texp (float): The exponent for the Box-Cox transformation. Default is 4.
    
    Returns:
    - ndarray: Transformed data.
    """
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    data[data<=0] = 0
    datat = (data ** (1 / texp) - 1) / (1 / texp)

    return datat

def boxcox_back_transform(data, texp=4):
    """
    Apply the inverse of the Box-Cox transformation to the input data.
    
    Parameters:
    - data (array-like): Transformed data to be reverted.
    - texp (float): The exponent used in the original Box-Cox transformation. Default is 4.
    
    Returns:
    - ndarray: Original (reverted) data.
    """
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    data[data<-texp] = -texp
    datat = (data / texp + 1) ** texp
    return datat

def lognormal_transform(data, c=0.001):
    """
    Apply the lognormal transformation to the input data.
    
    Parameters:
    - data (array-like): Input data to be transformed.
    - c (float): A small constant added to handle zero values. Default is 0.001.
    
    Returns:
    - ndarray: Transformed data.
    """
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    data[data < 0] = 0
    datat = np.log(data + c)
    return datat

def lognormal_back_transform(data, c=0.001):
    """
    Apply the inverse of the lognormal transformation to the input data.
    
    Parameters:
    - data (array-like): Transformed data to be reverted.
    - c (float): The constant used in the original lognormal transformation. Default is 0.001.
    
    Returns:
    - ndarray: Original (reverted) data.
    """
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    datat = np.exp(data) - c
    return datat

def empirical_cdf_transform(data,times,cdfs):
    """
    Apply the normal quantile transform to the input precipitation data for each month.

    This function transforms the empirical distribution of monthly precipitation data 
    to a standard normal distribution using the cumulative distribution function (CDF).

    Parameters:
    - data (array-like): 2D array of precipitation data with shape (nstn, ntime).
    - time (array-like): Array of dates corresponding to the precipitation data.

    Returns:
    - ndarray: Transformed data with Z-scores, shape (nstn, ntime).
    """
    if data.shape[0] != len(times):
        data = data.T

    df = pd.DataFrame(data=data, index=pd.to_datetime(times))
    transformed_df = normal_quantile_transform(df, cdfs)
    datat = transformed_df.to_numpy()

    datat = datat.T

    return datat
 
def empirical_cdf_back_transform(data, times, cdfs):
    """
    Apply the inverse of the normal quantile transform to the input transformed data.

    This function reverts the transformed Z-scores back to the original precipitation values
    using the inverse of the empirical cumulative distribution function (CDF).

    Parameters:
    - data (array-like): 2D array of transformed data (Z-scores) with shape (nstn, ntime).
    - time (array-like): Array of dates corresponding to the transformed data.

    Returns:
    - ndarray: Original (reverted) precipitation data, shape (nstn, ntime).
    """
    if data.shape[0] != len(times):
        data = data.T

    df = pd.DataFrame(data=data, index=pd.to_datetime(times))
    back_transformed_df = inverse_normal_quantile_transform(df, cdfs)
    datat = back_transformed_df.to_numpy()
    #datat = datat.T

    return datat


def data_transformation(data,method, settings, mode='transform',times=None,cdfs=None):
    """
    Transform or back-transform the input data based on the specified method and mode.
    
    Parameters:
    - data (array-like): Input data to be processed.
    - method (str): The transformation method. Accepts 'boxcox','lognormal' or 'monthly_cdf'.
    - settings (dict): Dictionary containing parameters for the transformation methods.
    - mode (str): Specifies the operation mode. Accepts 'transform' or 'back_transform'. Default is 'transform'.
    
    Returns:
    - ndarray: Processed data.
    """

    if method == 'boxcox':
        if mode == 'transform':
            data = boxcox_transform(data, settings['exponent'])
        elif mode == 'back_transform':
            data = boxcox_back_transform(data, settings['exponent'])
        else:
            print('Unknown transformation mode: entry=', mode)
            sys.exit()

    elif method == 'lognormal':
        if mode == 'transform':
            data = lognormal_transform(data, settings.get('constant', 0.001))
        elif mode == 'back_transform':
            data = lognormal_back_transform(data, settings.get('constant', 0.001))
        else:
            print('Unknown transformation mode: entry=', mode)
            sys.exit()
    
    elif method == 'empirical_cdf':
        if mode == 'transform':
            data = empirical_cdf_transform(data,times,cdfs)
        elif mode == 'back_transform':
            data = empirical_cdf_back_transform(data,times,cdfs)
        else:
            print('Unknown transformation mode: entry=', mode)
            sys.exit()
    else:
        print('Unknown transformation method: entry=', method)
        sys.exit()
    
    return data



########################################################################################################################
# input station data processing

def merge_stndata_into_single_file(config):
    # GMET v2.0 assumes that each station has an independent file. GPEP will merge the station data into one file,
    #   which can speed up i/o in subsequent runs of GPEP using the same dataset. 
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
        maxRange_vars = [np.inf] * len(target_vars)

    if 'transform_vars' in config:
        transform_vars = config['transform_vars']
        if not isinstance(transform_vars, list):
            transform_vars = [transform_vars] * len(target_vars)
    else:
        transform_vars = [''] * len(target_vars)

    if 'transform' in config:
        transform_settings = config['transform']
    else:
        transform_settings = {}

    if 'mapping_InOut_var' in config:
        mapping_InOut_var  = config['mapping_InOut_var']
    else:
        mapping_InOut_var = []

        
    if 'overwrite_merged_stnfile' in config:
        overwrite_merged_stnfile = config['overwrite_merged_stnfile']
    else:
        overwrite_merged_stnfile = True

    # settings and prints
    print('#' * 50)
    print('Merging individual station files to one single file')
    print('#' * 50)
    print('Input station list:     ', input_stn_list)
    print('Input station folder:   ', input_stn_path)
    print('Output station file:    ', file_allstn)
    print('Target variables:       ', input_vars)

    if os.path.isfile(file_allstn):
        print('NOTE: Merged station file exists')
        if overwrite_merged_stnfile == True:
            print('overwrite_merged_stnfile is True. Continue.')
        else:
            print('overwrite_merged_stnfile is False. Skip station merging.\n')
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
            if transform_vars[i] == 'empirical_cdf':
                cdfs = calculate_monthly_cdfs(config['file_allstn'], target_vars[i])
                ds_stn[tvar].values = data_transformation(ds_stn[target_vars[i]].values, transform_vars[i],
                                                transform_settings[transform_vars[i]], 'transform',
                                                times=ds_stn['time'].values,cdfs=cdfs)
            else:
                ds_stn[tvar].values = data_transformation(ds_stn[target_vars[i]].values, transform_vars[i],
                                                transform_settings[transform_vars[i]], 'transform')
        else:
            print(f'Do not perform transformation for {target_vars[i]}')

    # save to output files
    encoding = {}
    for var in ds_stn.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}

    if 'stnid' in ds_stn:
        try:
            ds_stn['stnid'] = ds_stn['stnid'].astype('|S')
        except:
            print('Failed to process station ID information.')

    ds_stn.to_netcdf(file_allstn, encoding=encoding)

    t2 = time.time()
    print('Time cost (s) for merging station data file:', t2-t1, '\n')

    return config


