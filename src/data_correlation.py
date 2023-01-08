# calculate temporal auto correlation or spatial correlation length
import numpy as np
import xarray as xr
from near_stn_search import distance
from scipy.optimize import curve_fit
import os


def cal_cross_cc(d1, d2):
    index = ~np.isnan(d1) & ~np.isnan(d2)
    if np.sum(index) < 2:
        cc = np.nan
    else:
        cc = np.corrcoef(d1[index], d2[index])[0, 1]
    return cc


def cal_auto_cc(data, lag=0):
    cc = cal_cross_cc(data[:-lag], data[lag:])
    return cc


def station_lag_correlation(data, lag):
    # data: [number of stations, number of time steps]
    nstn, ntime = data.shape
    cc = np.nan * np.zeros(nstn)
    for i in range(nstn):
        cc[i] = cal_auto_cc(data[i, :], lag)
    return cc


def station_cross_correlation(data1, data2):
    # data1/data2: [number of stations, number of time steps]
    nstn, ntime = data1.shape
    cc = np.nan * np.zeros(nstn)
    for i in range(nstn):
        cc[i] = cal_cross_cc(data1[i, :], data2[i, :])
    return cc


def station_spatial_correlation(lat, lon, data):
    # cc between each pair of stations
    nstn = len(lat)
    npair = int(((nstn - 1) * nstn / 2))
    cc_pair = np.nan * np.zeros(npair)
    dist_pair = np.nan * np.zeros(npair)
    index_pair = np.nan * np.zeros([npair, 2])

    flag = 0
    for i in range(nstn-1):
        for j in range(i+1, nstn):
            cc_pair[flag] = cal_cross_cc(data[i, :], data[j, :])
            dist_pair[flag] = distance(lat[i], lon[i], lat[j], lon[j])
            index_pair[flag, 0] = i
            index_pair[flag, 1] = j
            flag = flag + 1

    return cc_pair, dist_pair, index_pair



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
    # parse and change configurations
    path_stn_info = config['path_stn_info']
    file_stn_cc = f'{path_stn_info}/stn_time_space_correlation.nc'
    config['file_stn_cc'] = file_stn_cc

    # in/out information to this function
    file_allstn = config['file_allstn']
    target_vars = config['target_vars']
    file_stn_cc = config['file_stn_cc']
    overwrite_station_cc = config['overwrite_station_cc']

    lag = 1 # lag-1 cc

    print('#' * 50)
    print(f'Calculate station space and time correlation')
    print('#' * 50)
    print('Input file_allstn:', file_allstn)
    print('Output file_stn_cc:', file_stn_cc)
    print('Target variables:', target_vars)

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
            stn_value = ds_stn[var_name].values
            lat = ds_stn['lat'].values
            lon = ds_stn['lon'].values

        cc = station_lag_correlation(stn_value, lag)

        ds_out[var_name + '_cc_lag1'] = xr.DataArray(cc, dims=('stn'))
        ds_out[var_name + '_cc_lag1_mean'] = xr.DataArray([np.nanmean(cc)], dims=('z'))

        # calculate space correlation
        print('Spatial correlation calculation for:', var_name)
        cc_pair, dist_pair, index_pair = station_spatial_correlation(lat, lon, stn_value)

        index = ~np.isnan(dist_pair + cc_pair)
        # popt, pcov = curve_fit(func_clen_exp2p, dist_pair[index], cc_pair[index])
        popt, pcov = curve_fit(func_clen_exp1p, dist_pair[index], cc_pair[index])

        # ds_out[f'{var_name}_space_pair_cc'] = xr.DataArray(cc_pair, dims=('pair'))
        # ds_out[f'{var_name}_space_pair_dist'] = xr.DataArray(dist_pair, dims=('pair'))
        # ds_out[f'{var_name}_space_pair_index'] = xr.DataArray(index_pair, dims=('pair', 'pind'))
        ds_out[f'{var_name}_space_Clen'] = xr.DataArray([popt[0]], dims=('z'))
        # ds_out[f'{var_name}_space_param2'] = xr.DataArray([popt[1]], dims=('z'))

    ########################################################################################################################
    # calculate cross correlation
    # just precipitation VS trange. This should be changed in the future to variable dependence setting without any hard-caoded variable

    if ('prcp' in target_vars) and ('trange' in target_vars):
        print('Calculate cross correlation between prcp and trange')
        with xr.open_dataset(file_allstn) as ds_stn:
            data1 = ds_stn['prcp'].values
            data2 = ds_stn['trange'].values

        cc = station_cross_correlation(data1, data2)

        ds_out['prcp_trange_cc_cross'] = xr.DataArray(cc, dims=('stn'))
        ds_out['prcp_trange_cc_cross_mean'] = xr.DataArray([np.nanmean(cc)], dims=('z'))

    ########################################################################################################################
    # save output file
    encoding = {}
    for var in ds_out.data_vars:
        encoding[var] = {'zlib': True, 'complevel': 4}

    ds_out.to_netcdf(file_stn_cc, encoding=encoding)

    return config