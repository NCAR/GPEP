# calculate temporal auto correlation or spatial correlation length
import sys, time

import numpy as np
import xarray as xr
from near_stn_search import distance
from scipy.optimize import curve_fit
from multiprocessing import Pool
import os


def cal_cross_cc(d1, d2):
    index = ~np.isnan(d1) & ~np.isnan(d2)
    if np.sum(index) < 2:
        cc = np.nan
    else:
        cc = np.corrcoef(d1[index], d2[index])[0, 1]
    return cc


def station_lag_correlation(data, lag):
    # data: [number of stations, number of time steps]
    nstn, ntime = data.shape
    cc = np.nan * np.zeros(nstn)
    for i in range(nstn):
        d1 = data[i, :-lag]
        d2 = data[i, lag:]
        ind = ~np.isnan(d1+d2)
        if np.sum(ind) > 2:
            cc[i] = cal_cross_cc(d1[ind], d2[ind])
    return cc


def station_cross_correlation(data1, data2):
    # data1/data2: [number of stations, number of time steps]
    nstn, ntime = data1.shape
    cc = np.nan * np.zeros(nstn)
    for i in range(nstn):
        cc[i] = cal_cross_cc(data1[i, :], data2[i, :])
    return cc

# parallel version
def init_worker(data, lat, lon):
    # Using a dictionary is not strictly necessary. You can also
    # use global variables.
    global mppool_ini_dict
    mppool_ini_dict = {}
    mppool_ini_dict['data'] = data
    mppool_ini_dict['lat'] = lat
    mppool_ini_dict['lon'] = lon

def single_cc_dist(i, j):
    cc_pair = cal_cross_cc(mppool_ini_dict['data'][i, :], mppool_ini_dict['data'][j, :])
    dist_pair = distance(mppool_ini_dict['lat'][i], mppool_ini_dict['lon'][i], mppool_ini_dict['lat'][j], mppool_ini_dict['lon'][j])
    return cc_pair, dist_pair


def station_spatial_correlation(lat, lon, data, num_processes=4):
    # cc between each pair of stations
    nstn = len(lat)
    npair = int(((nstn - 1) * nstn / 2))
    cc_pair = np.nan * np.zeros(npair)
    dist_pair = np.nan * np.zeros(npair)
    index_pair = np.nan * np.zeros([npair, 2])


    items = [(i, j) for i in range(nstn - 1) for j in range(i + 1, nstn)]
    with Pool(processes=num_processes, initializer=init_worker, initargs=(data, lat, lon)) as pool:
        result = pool.starmap(single_cc_dist, items)

    flag = 0
    for i in range(len(items)):
        cc_pair[flag] = result[i][0]
        dist_pair[flag] = result[i][1]
        index_pair[flag, 0] = items[i][0]
        index_pair[flag, 1] = items[i][1]
        flag = flag + 1

    return cc_pair, dist_pair, index_pair

# # serial version
# def station_spatial_correlation(lat, lon, data):
#     # cc between each pair of stations
#     nstn = len(lat)
#     npair = int(((nstn - 1) * nstn / 2))
#     cc_pair = np.nan * np.zeros(npair)
#     dist_pair = np.nan * np.zeros(npair)
#     index_pair = np.nan * np.zeros([npair, 2])
#
#     flag = 0
#     for i in range(nstn-1):
#         for j in range(i+1, nstn):
#             cc_pair[flag] = cal_cross_cc(data[i, :], data[j, :])
#             dist_pair[flag] = distance(lat[i], lon[i], lat[j], lon[j])
#             index_pair[flag, 0] = i
#             index_pair[flag, 1] = j
#             flag = flag + 1
#
#     return cc_pair, dist_pair, index_pair



def func_clen_exp1p(x, a):
    # regression of spatial correlation length
    y = np.exp(-x/a)
    return y

def func_clen_Gaussian_gstools(x, a):
    # regression of spatial correlation length
    y = np.exp(-(x/a)**2)
    return y


# def func_clen_exp2p(x, a, b):
#     # regression of spatial correlation length
#     y = np.exp(-(x/a)**b)
#     return y


def func_clen_exp2p(x, a, b):
    # regression of spatial correlation length
    y = b * np.exp(-x/a)
    return y


def func_clen_exp3p(x, a, b, c):
    # regression of spatial correlation length
    y = c * np.exp(- (x/a) ** b)
    return y


def hist_xy(x, y):
    binnum = 20 # divide x into binnum bins
    smpnum = int(len(x)/binnum) # sample number in each bin
    index = np.argsort(x)
    x = x[index]
    y = y[index]
    x2 = []
    y2 = []

    for i in range(binnum):
        x2.append(np.median(x[i * smpnum:(i + 1) * smpnum]))
        y2.append(np.median(y[i * smpnum:(i + 1) * smpnum]))

    return np.array(x2), np.array(y2)


def station_space_time_correlation(config):
    t1 = time.time()

    # parse and change configurations
    path_stn_info = config['path_stn_info']
    file_stn_cc = f'{path_stn_info}/stn_time_space_correlation.nc'
    config['file_stn_cc'] = file_stn_cc

    # in/out information to this function
    file_allstn = config['file_allstn']
    target_vars = config['target_vars']
    file_stn_cc = config['file_stn_cc']
    num_processes = config['num_processes']
    auto_corr_method = config['auto_corr_method']
    rolling_window = config['rolling_window']

    if 'clen' in config:
        clen_config = config['clen']
        if not isinstance(clen_config, list):
            clen_config = [clen_config] * len(target_vars)
    else:
        clen_config = [-9999] * len(target_vars)

    if 'overwrite_station_cc' in config:
        overwrite_station_cc = config['overwrite_station_cc']
    else:
        overwrite_station_cc = False

    linkvar0 = config['linkvar']
    linkvar = {}
    for v in linkvar0:
        linkvar[v[0]] = v[1]

    lag = 1 # lag-1 cc

    print('#' * 50)
    print(f'Calculate station space and time correlation')
    print('#' * 50)
    print('Input file_allstn:', file_allstn)
    print('Output file_stn_cc:', file_stn_cc)
    print('Target variables:', target_vars)
    print('Number of processes:', num_processes)
    print('Linked variables:', linkvar)

    if os.path.isfile(file_stn_cc):
        print('Note! Station correlation file exists')
        if overwrite_station_cc == True:
            print('overwrite_station_cc is True. Continue.')
        else:
            print('overwrite_station_cc is False. Skip correlation calculation.')
            return config

    ########################################################################################################################
    # initialize outputs
    with xr.open_dataset(file_allstn) as ds_out:
        for v in ds_out.data_vars:
            if 'time' in ds_out[v].dims:
                ds_out = ds_out.drop_vars(v)

    ds_out.coords['z'] = np.arange(1)
    ds_out.coords['pair'] = np.arange(int(len(ds_out.stn) * (len(ds_out.stn)-1) / 2))
    ds_out.coords['pind'] = np.arange(2)

    ########################################################################################################################

    for vn in range(len(target_vars)):

        var_name = target_vars[vn]

        # calculate auto correlation
        print('Auto correlation calculation for:', var_name)
        with xr.open_dataset(file_allstn) as ds_stn:

            # calculate var - moving_averaging(var) to remove monthly cycle
            if auto_corr_method == 'direct':
                stn_value = ds_stn[var_name].values
            elif auto_corr_method == 'anomaly':
                stn_value_raw = ds_stn[var_name].values
                stn_value_mv = ds_stn[var_name].rolling(time=rolling_window, center=True).mean().values
                stn_value = stn_value_raw - stn_value_mv
            else:
                sys.exit('Unknown auto_corr_method')

            lat = ds_stn['lat'].values
            lon = ds_stn['lon'].values

        cc = station_lag_correlation(stn_value, lag)

        ds_out[var_name + '_cc_lag1'] = xr.DataArray(cc, dims=('stn'))
        ds_out[var_name + '_cc_lag1_mean'] = xr.DataArray([np.nanmean(cc)], dims=('z'))

        # calculate space correlation
        if not clen_config[vn] > 0:
            print('Spatial correlation calculation for:', var_name)
            cc_pair, dist_pair, index_pair = station_spatial_correlation(lat, lon, stn_value, num_processes)

            index = ~np.isnan(dist_pair + cc_pair)
            # popt, pcov = curve_fit(func_clen_exp2p, dist_pair[index], cc_pair[index])
            popt, pcov = curve_fit(func_clen_exp1p, dist_pair[index], cc_pair[index])

            ds_out[f'{var_name}_space_pair_cc'] = xr.DataArray(cc_pair, dims=('pair'))
            ds_out[f'{var_name}_space_pair_dist'] = xr.DataArray(dist_pair, dims=('pair'))
            ds_out[f'{var_name}_space_pair_index'] = xr.DataArray(index_pair, dims=('pair', 'pind'))
            ds_out[f'{var_name}_space_Clen'] = xr.DataArray([popt[0]], dims=('z'))
            # ds_out[f'{var_name}_space_param2'] = xr.DataArray([popt[1]], dims=('z'))

    ########################################################################################################################
    # calculate cross correlation
    # just precipitation VS trange. This should be changed in the future to variable dependence setting without any hard-caoded variable

    # calculate var - moving_averaging(var) to remove monthly cycle
    for var1, var2 in linkvar.items():
        print(f'Calculate cross correlation between {var1} and {var2}')
        with xr.open_dataset(file_allstn) as ds_stn:
            # data1 = ds_stn['prcp'].values

            if auto_corr_method == 'direct':
                data1 = ds_stn[var1].values
                data2 = ds_stn[var2].values
            elif auto_corr_method == 'anomaly':
                stn_value_raw = ds_stn[var1].values
                stn_value_mv = ds_stn[var1].rolling(time=rolling_window, center=True).mean().values
                data1 = stn_value_raw - stn_value_mv

                stn_value_raw = ds_stn[var2].values
                stn_value_mv = ds_stn[var2].rolling(time=rolling_window, center=True).mean().values
                data2 = stn_value_raw - stn_value_mv
            else:
                sys.exit('Unknown auto_corr_method')

        cc = station_cross_correlation(data1, data2)

        ds_out[f'{var1}_{var2}_cc_cross'] = xr.DataArray(cc, dims=('stn'))
        ds_out[f'{var1}_{var2}_cc_cross_mean'] = xr.DataArray([np.nanmean(cc)], dims=('z'))

    ########################################################################################################################
    # save output file
    encoding = {}
    for var in ds_out.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}

    ds_out.to_netcdf(file_stn_cc, encoding=encoding)

    t2 = time.time()
    print('Time cost (seconds):', t2 - t1)
    print('Successful correlation calculation!\n\n')

    return config