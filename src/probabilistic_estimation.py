### routines to generate ensembles (probabilistic estimation) from regression outputs

# load libraries
import os, sys, time, numbers
import xarray as xr
import numpy as np
from scipy import special
from multiprocessing import Pool

import random_field_FortranGMET as rf_FGMET
from data_processing import data_transformation
from empirical_cdf import calculate_monthly_cdfs

# ====== subroutines/ functions ======

def perturb_estimates_withoccurrence(data, uncert, poe, rndnum, minrndnum=-3.99, maxrndnum=3.99):
    data     = data.copy()
    uncert   = uncert.copy()
    poe      = poe.copy()
    rndnum   = rndnum.copy()
    
    # initialize output array
    data_out = np.nan * np.zeros(data.shape, dtype=np.float32)

    #### calculate conditional precipitation probability: cprob
    acorr = rndnum / np.sqrt(2)
    aprob = special.erfc(acorr)
    cprob = (2 - aprob) / 2
    cprob[cprob < 3e-5] = 3e-5
    cprob[cprob > 0.99997] = 0.

    ##### positive precipitation
    index_positive = cprob >= (1 - poe)
    cs = (cprob - (1 - poe)) / poe
    cs[~index_positive] = np.nan  # because poe == 0 is not excluded in the matrix
    
    # assign a small precipitation value to positive grids if their values = 0
    dtemp = data[index_positive]
    dtemp[dtemp < 0.1] = 0.1
    data[index_positive] = dtemp
    del dtemp
    
    # generate random numbers
    cs[cs > 0.99997] = 0.99997
    cs[cs < 3e-5] = 3e-5
    rn = np.sqrt(2) * special.erfinv(2 * cs - 1)
    rn[rn < minrndnum] = minrndnum
    rn[rn > maxrndnum] = maxrndnum

    ## start probabilistic estimation
    # generate mu and sigma from lognormal distribution
    data_out[index_positive] = data[index_positive] + rn[index_positive] * uncert[index_positive]

    ## zero precipitation
    index_nonpositive = cprob < (1 - poe)
    data_out[index_nonpositive] = np.nanmin(data_out)
    del index_nonpositive, index_positive

    return data_out


def perturb_estimates_general(data, uncert, rndnum, minrndnum=-3.99, maxrndnum=3.99):
    data   = data.copy()
    uncert = uncert.copy()
    rndnum = rndnum.copy()

    rndnum[rndnum < minrndnum] = minrndnum
    rndnum[rndnum > maxrndnum] = maxrndnum

    data_out = data + rndnum * uncert
    return data_out


def prob_estimate_for_one_var(var_name, reg_estimate, reg_error, nearby_stn_max, poe, random_field, minrndnum, maxrndnum, 
                              transform_method, transform_setting, ds_out,config):
    # generate probabilistic estimates
    if poe.shape == reg_estimate.shape:
        ens_estimate = perturb_estimates_withoccurrence(reg_estimate, reg_error, poe, random_field, minrndnum, maxrndnum)
        # back transformation
        if len(transform_method) > 0:
            if transform_method == 'empirical_cdf':
                cdfs = calculate_monthly_cdfs(xr.open_dataset(config['file_allstn']),var_name)
                ens_estimate = data_transformation(ens_estimate, transform_method, transform_setting, 'back_transform',
                                                times=ds_out['time'].values,cdfs=cdfs)
            else:
                ens_estimate = data_transformation(ens_estimate, transform_method, transform_setting, 'back_transform')
            # variable and unit dependent ...
            # minprcp = 0.01
            # ens_estimate[ens_estimate < minprcp] = 0
            # ens_estimate[(ens_estimate > minprcp) & (ens_estimate < minprcp)] = minprcp
    else:
        ens_estimate = perturb_estimates_general(reg_estimate, reg_error, random_field, minrndnum, maxrndnum)

    # max limit: may be changed in the future ...
    if np.array(nearby_stn_max).shape == reg_estimate.shape:
        if len(transform_method) > 0:
            precip_err_cap = 0.2 # hard coded ...
            if transform_method == 'empirical_cdf':
                cdfs = calculate_monthly_cdfs(xr.open_dataset(config['file_allstn']),var_name)
                nearby_stn_max = data_transformation(nearby_stn_max+reg_error*precip_err_cap, transform_method, transform_setting, 'back_transform',
                                                times=ds_out['time'].values,cdfs=cdfs)
            else:
                nearby_stn_max = data_transformation(nearby_stn_max+reg_error*precip_err_cap, transform_method, transform_setting, 'back_transform')
        else:
            nearby_stn_max = nearby_stn_max + reg_error * 2
        mask = ens_estimate > nearby_stn_max
        ens_estimate[mask] = nearby_stn_max[mask]

    # # add to output ds
    ds_out[var_name] = xr.DataArray(ens_estimate, dims=('y', 'x', 'time'))

    return ds_out


def generate_random_numbers(masterseed, n=1):
    np.random.seed(masterseed)
    return np.random.randint(0, 2**16, n)


def spcorr_structure(config):

    t1 = time.time()

    # parse and change configurations
    outpath_parent            = config['outpath_parent']
    path_ensemble             = f'{outpath_parent}/ensembles/'
    os.makedirs(path_ensemble, exist_ok=True)
    config['path_ensemble']   = path_ensemble
    file_ens_prefix           = f'{path_ensemble}/Ensemble_estimate_'    # this may get overwritten before use
    config['file_ens_prefix'] = file_ens_prefix


    path_spcorr = f'{outpath_parent}/spcorr'
    os.makedirs(path_spcorr, exist_ok=True)
    config['path_rfweight'] = path_spcorr

    file_spcorr_prefix = f'{path_spcorr}/spcorr_'
    config['file_spcorr_prefix'] = file_spcorr_prefix

    file_stn_nearinfo = config['file_stn_nearinfo']

    # in/out information to this function
    file_stn_cc        = config['file_stn_cc']
    file_grid_reg      = config['file_grid_reg']
    file_spcorr_prefix = file_spcorr_prefix
    target_vars        = config['target_vars']
    grid_lat_name      = config['grid_lat_name']
    grid_lon_name      = config['grid_lon_name']

    if 'clen' in config:
        clen_config = config['clen']
        if not isinstance(clen_config, list):
            clen_config = [clen_config] * len(target_vars)
    else:
        clen_config = [-9999] * len(target_vars)

    if 'overwrite_spcorr' in config:
        overwrite_spcorr = config['overwrite_spcorr']
    else:
        overwrite_spcorr = False

    print(f'Generating spatial correlation structure')

    # spatial correlation
    allvar_clen = {}
    with xr.open_dataset(file_stn_cc) as ds_stn_cc:
        for i in range(len(target_vars)):
            var_name = target_vars[i]
            if clen_config[i] > 0:
                allvar_clen[var_name] = clen_config[i]
            else:
                allvar_clen[var_name] = ds_stn_cc[var_name + '_space_Clen'].values

    # regression estimates
    with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
        grid_lat = ds_nearinfo[grid_lat_name].values
        grid_lon = ds_nearinfo[grid_lon_name].values

    for var_name in target_vars:
        file_spcor = f'{file_spcorr_prefix}{var_name}.npz'
        if os.path.isfile(file_spcor):
            print(f'spcorr outfile exists: {file_spcor}')
            if overwrite_spcorr == True:
                print('overwrite_spcorr is True. Overwrite it.')
                _ = os.remove(file_spcor)
                rf_FGMET.spcorr_grd(grid_lat, grid_lon, allvar_clen[var_name], file_spcor)
            else:
                print('overwrite_spcorr is False. Skip spcorr generation.\n')
                continue
        else:
            print('Creating spcorr outfile:', file_spcor)
            rf_FGMET.spcorr_grd(grid_lat, grid_lon, allvar_clen[var_name], file_spcor)

    t2 = time.time()
    print('Successful structure generation. Time cost (sec):', t2-t1)


def generate_prob_estimates_serial(config, member_range=[]):
    # most time cost comes from field_rand

    t1 = time.time()

    # parse and change configurations
    case_name                 = config['case_name']
    outpath_parent            = config['outpath_parent']
    path_ensemble             = f'{outpath_parent}/ensembles/'
    os.makedirs(path_ensemble, exist_ok=True)
    config['path_ensemble']   = path_ensemble
    file_ens_prefix           = f'{path_ensemble}/{case_name}_ensMember_'
    config['file_ens_prefix'] = file_ens_prefix

    path_spcorr                  = f'{outpath_parent}/spcorr'
    os.makedirs(path_spcorr, exist_ok=True)
    config['path_rfweight']      = path_spcorr
    file_spcorr_prefix           = f'{path_spcorr}/spcorr_'
    config['file_spcorr_prefix'] = file_spcorr_prefix

    # in/out information to this function
    file_stn_cc         = config['file_stn_cc']
    file_grid_reg       = config['file_grid_reg']
    file_grid_auxiliary = config['file_grid_auxiliary']
    #file_ens_prefix     = config['file_ens_prefix']   # redundant
    #file_spcorr_prefix  = file_spcorr_prefix
    ensemble_start      = config['ensemble_start']
    ensemble_end        = config['ensemble_end']
    target_vars         = config['target_vars']

    if 'target_vars_WithProbability' in config:
        target_vars_WithProbability = config['target_vars_WithProbability']
    else:
        target_vars_WithProbability = []

    if 'probability_thresholds' in config:
        probability_thresholds      = config['probability_thresholds']
    else:
        probability_thresholds = [0] * len(target_vars_WithProbability)

    file_stn_nearinfo = config['file_stn_nearinfo']

    if 'master_seed' in config:
        master_seed = config['master_seed']
    else:
        master_seed = -1
    
    if 'append_date_to_output_filename' in config:
        append_date_to_output_filename = config['append_date_to_output_filename']
    else:
        append_date_to_output_filename = True
        
    datestamp = f"{config['date_start'].replace('-', '')}-{config['date_end'].replace('-', '')}"

    grid_lat_name = config['grid_lat_name']
    grid_lon_name = config['grid_lon_name']

    if 'transform_vars' in config:
        transform_vars = config['transform_vars']
    else:
        transform_vars = [''] * len(target_vars)

    if 'transform' in config:
        transform_settings = config['transform']
    else:
        transform_settings = {}

    if 'target_vars_max_constrain' in config:
        target_vars_max_constrain = config['target_vars_max_constrain']
    else:
        target_vars_max_constrain = []

    if 'output_randomfield' in config:
        output_randomfield = config['output_randomfield']
    else:
        output_randomfield = False

    if 'overwrite_ens' in config:
        overwrite_ens = config['overwrite_ens']
    else:
        overwrite_ens = False

    linkvar0 = config['linkvar']

    if 'clen' in config:
        clen_config = config['clen']
        if not isinstance(clen_config, list):
            clen_config = [clen_config] * len(target_vars)
    else:
        clen_config = [-9999] * len(target_vars)

    if 'lag1_auto_cc' in config:
        lag1_auto_cc = config['lag1_auto_cc']
        if not isinstance(lag1_auto_cc, list):
            lag1_auto_cc = [lag1_auto_cc] * len(target_vars)
    else:
        lag1_auto_cc = [-9999] * len(target_vars)

    if 'cross_cc' in config:
        cross_cc = config['cross_cc']
        if not isinstance(cross_cc, list):
            cross_cc = [cross_cc] * len(linkvar0)
    else:
        cross_cc = [-9999] * len(linkvar0)

    maxrndnum = 3.99  # minimum and maximum random number following GMET scripts
    minrndnum = -3.99

    linkvar = {}
    for i in range(len(linkvar0)):
        v = linkvar0[i]
        linkvar[v[0]] = v[1]
        linkvar[f"{v[0]}_crosscc"] = cross_cc[i]

    ensemble_number = ensemble_end - ensemble_start + 1
    if len(member_range) == 0:
        member_range = [0, ensemble_number]

    member_range_2 = np.array(member_range) + ensemble_start
    print(f'Starting probabilistic estimation of {target_vars} with linked vars {linkvar} for member range {member_range_2}')

    ########################################################################################################################
    # create all seeds used to generate random fields
    with xr.open_dataset(file_grid_reg) as ds_grid_reg:
        date0 = int(ds_grid_reg.time[0].dt.strftime('%Y%m%d%H'))
        ntime = len(ds_grid_reg.time)

    if master_seed >= 0:
        master_seed = master_seed + date0 # ensure different input batches have different seeds
    else:
        master_seed = np.random.randint(1e9)

    seeds_rf  = generate_random_numbers(master_seed, len(target_vars) * ensemble_number * ntime)
    seeds_rf  = np.reshape(seeds_rf, [len(target_vars), ensemble_number, ntime])
    seeds_rf2 = {}
    for i in range(len(target_vars)):
        seeds_rf2[target_vars[i]] = seeds_rf[i, :, :]
    seeds_rf = seeds_rf2
    del seeds_rf2

    ########################################################################################################################
    # load data for regression

    allvar_auto_lag1_cc = {}
    allvar_clen         = {}
    allvar_reg_estimate = {}
    allvar_reg_error    = {}
    allvar_poe          = {}
    nearby_stn_max      = {}

    # auto correlation
    with xr.open_dataset(file_stn_cc) as ds_stn_cc:
        for i in range(len(target_vars)):
            var_name = target_vars[i]

            if lag1_auto_cc[i] >= -1:
                allvar_auto_lag1_cc[var_name] = lag1_auto_cc[i]
            else:
                allvar_auto_lag1_cc[var_name] = ds_stn_cc[var_name + '_cc_lag1_mean'].values

            if clen_config[i] > 0:
                allvar_clen[var_name] = clen_config[i]
            else:
                allvar_clen[var_name] = ds_stn_cc[var_name + '_space_Clen'].values

    # regression estimates
    with xr.open_dataset(file_grid_reg) as ds_grid_reg:
        for vn in range(len(target_vars)):
            var_name = target_vars[vn]
            if len(transform_vars[vn]) > 0:
                allvar_reg_estimate[var_name] = ds_grid_reg[var_name + '_' + transform_vars[vn]].values
            else:
                allvar_reg_estimate[var_name] = ds_grid_reg[var_name].values

            if var_name in target_vars_WithProbability:
                var_poe = var_name + '_poe'
                allvar_poe[var_name] = ds_grid_reg[var_poe].values
            else:
                allvar_poe[var_name] = np.array([])

        tartime = ds_grid_reg.time.values

    # auxiliary info: estimate error and nearby_stn_max
    with xr.open_dataset(file_grid_auxiliary) as ds_grid_aux:
        for vn in range(len(target_vars)):
            var_name = target_vars[vn]
            if len(transform_vars[vn]) > 0:
                allvar_reg_error[var_name] = ds_grid_aux['uncert_' + var_name + '_' + transform_vars[vn]].values
                if var_name in target_vars_max_constrain:
                    nearby_stn_max[var_name] = ds_grid_aux['nearmax_' + var_name + '_' + transform_vars[vn]].values
                else:
                    nearby_stn_max[var_name] = np.array([])
            else:
                allvar_reg_error[var_name] = ds_grid_aux['uncert_' + var_name].values
                if var_name in target_vars_max_constrain:
                    nearby_stn_max[var_name] = ds_grid_aux['nearmax_' + var_name].values
                else:
                    nearby_stn_max[var_name] = np.array([])

    nrow, ncol, ntime = allvar_reg_error[var_name].shape

    ########################################################################################################################
    # dependent/independent variables and their correlation

    target_vars_independent = []
    target_vars_dependent = []
    target_vars_dependent_cross_cc = []
    for v in target_vars:
        if not v in linkvar:
            target_vars_independent.append(v)
        else:
            target_vars_dependent.append(v)

            if linkvar[f'{v}_crosscc'] >= -1:
                crosscc = linkvar[f'{v}_crosscc']
            else:
                cross_cc_varname1 = f'{v}_{linkvar[v]}_cc_cross_mean'
                cross_cc_varname2 = f'{linkvar[v]}_{v}_cc_cross_mean'
                if cross_cc_varname1 in ds_stn_cc:
                    crosscc = ds_stn_cc[cross_cc_varname1].values[0]
                elif cross_cc_varname2 in ds_stn_cc:
                    crosscc = ds_stn_cc[cross_cc_varname1].values[0]
                else:
                    sys.exit(f'Cannot find {cross_cc_varname1} or {cross_cc_varname1} in {file_stn_cc}')

            target_vars_dependent_cross_cc.append(crosscc)

    # convert transform_vars to dict
    transform_vars2 = {}
    for i in range(len(target_vars)):
        transform_vars2[target_vars[i]] = transform_vars[i]
        if not transform_vars[i] in transform_settings:
            transform_settings[transform_vars[i]] = ''

    transform_vars = transform_vars2
    del transform_vars2

    ########################################################################################################################
    # generate spatial correlation structure

    spcorr_jpos = {}
    spcorr_ipos = {}
    spcorr_wght = {}
    spcorr_sdev = {}
    iorder      = {}
    jorder      = {}
    for var_name in target_vars:
        file_spcor = f'{file_spcorr_prefix}{var_name}.npz'
        dtmp = np.load(file_spcor, allow_pickle=True)
        spcorr_jpos[var_name] = dtmp['spcorr_jpos']
        spcorr_ipos[var_name] = dtmp['spcorr_ipos']
        spcorr_wght[var_name] = dtmp['spcorr_wght']
        spcorr_sdev[var_name] = dtmp['spcorr_sdev']
        iorder[var_name]      = dtmp['iorder']
        jorder[var_name]      = dtmp['jorder']

    ########################################################################################################################
    # generate ensemble members (probabilistic estimates)

    for ens in range(member_range[0], member_range[1]):
        # define output file name
        if append_date_to_output_filename == True:
            outfile_ens = f'{file_ens_prefix}{datestamp}_{ens + ensemble_start:03}.nc'
        else:
            outfile_ens = f'{file_ens_prefix}{ens + ensemble_start:03}.nc'
            
        if os.path.isfile(outfile_ens):
            print(f'Ensemble outfile exists: {outfile_ens}')
            if overwrite_ens == True:
                print('overwrite_ens is True. Overwrite it.')
            else:
                print('overwrite_ens is False. Skip probabilistic estimation.\n')
                continue

        # initialize outputs
        ds_out                = xr.Dataset()
        ds_out.coords['time'] = tartime

        with xr.open_dataset(file_stn_nearinfo) as ds_nearinfo:
            ds_out.coords['x']    = ds_nearinfo.coords['x']
            ds_out.coords['y']    = ds_nearinfo.coords['y']
            ds_out[grid_lat_name] = ds_nearinfo[grid_lat_name]
            ds_out[grid_lon_name] = ds_nearinfo[grid_lon_name]

        # loop over variables

        for vn in range(len(target_vars_independent)):

            var_name = target_vars_independent[vn]
            # print(f'Probabilistic estimation for {var_name} and ensemble member {ens}--{ensemble_number}')

            # generate random numbers
            random_field = np.nan * np.zeros([nrow, ncol, ntime], dtype=np.float32)
            for i in range(ntime):
                rndi = rf_FGMET.field_rand(spcorr_jpos[var_name], spcorr_ipos[var_name], spcorr_wght[var_name], spcorr_sdev[var_name],
                                           iorder[var_name], jorder[var_name], seeds_rf[var_name][ens, i])
                if i == 0:
                    random_field[:, :, i] = rndi
                else:
                    random_field[:, :, i] = random_field[:, :, i - 1] * allvar_auto_lag1_cc[var_name] + np.sqrt(1 - allvar_auto_lag1_cc[var_name] ** 2) * rndi

            random_field[np.isnan(allvar_reg_estimate[var_name])] = np.nan

            # probabilistic estimation (ensemble generation)
            ds_out = prob_estimate_for_one_var(var_name, allvar_reg_estimate[var_name], allvar_reg_error[var_name],
                                               nearby_stn_max[var_name], allvar_poe[var_name], random_field, 
                                               minrndnum, maxrndnum, transform_vars[var_name],
                                               transform_settings[transform_vars[var_name]], ds_out,config)
            
            if output_randomfield == True:
                ds_out[var_name + '_rnd'] = xr.DataArray(random_field, dims=('y', 'x', 'time'))

            # is any linked variable
            for d in range(len(target_vars_dependent)):
                var_name_dep = target_vars_dependent[d]
                if linkvar[var_name_dep] == var_name:
                    # print(f'Probabilistic estimation for {var_name_dep} and ensemble member {ens}--{ensemble_number}')

                    # generate random numbers
                    random_field_dep = np.nan * np.zeros([nrow, ncol, ntime], dtype=np.float32)
                    for i in range(ntime):
                        rndi = rf_FGMET.field_rand(spcorr_jpos[var_name_dep], spcorr_ipos[var_name_dep], spcorr_wght[var_name_dep],
                                                   spcorr_sdev[var_name_dep], iorder[var_name_dep], jorder[var_name_dep], 
                                                   seeds_rf[var_name_dep][ens, i])
                        random_field_dep[:, :, i] = rndi
                    random_field_dep = random_field * target_vars_dependent_cross_cc[d] + np.sqrt(1 - target_vars_dependent_cross_cc[d]**2) * random_field_dep
                    del random_field
                    random_field_dep[np.isnan(allvar_reg_estimate[var_name_dep])] = np.nan

                    # probabilistic estimation
                    ds_out = prob_estimate_for_one_var(var_name_dep, allvar_reg_estimate[var_name_dep], allvar_reg_error[var_name_dep],
                                                       nearby_stn_max[var_name], allvar_poe[var_name_dep], random_field_dep, 
                                                       minrndnum, maxrndnum, transform_vars[var_name_dep],
                                                       transform_settings[transform_vars[var_name_dep]], ds_out,config)

                    if output_randomfield == True:
                        ds_out[var_name + '_rnd'] = xr.DataArray(random_field_dep, dims=('y', 'x', 'time'))

        # save output file
        ds_out = ds_out.transpose('time', 'y', 'x')

        # reduce coordinate dims if needed
        if 'x' in ds_out.dims and 'y' in ds_out.dims:
            grid_lat_diff = np.abs(ds_out[grid_lat_name].isel(x=0).values - ds_out[grid_lat_name].isel(x=-1).values)
            grid_lon_diff = np.abs(ds_out[grid_lon_name].isel(y=0).values - ds_out[grid_lon_name].isel(y=-1).values)
            if (np.nanmax(grid_lat_diff) < 1e-10) and (np.nanmax(grid_lon_diff) < 1e-10):
                ds_out.coords['x'] = ds_out[grid_lon_name].isel(y=0).values
                ds_out.coords['y'] = ds_out[grid_lat_name].isel(x=0).values
                ds_out = ds_out.drop_vars([grid_lat_name, grid_lon_name])
                ds_out = ds_out.rename({'y': 'lat', 'x': 'lon'})

        ds_out = ds_out.fillna(-9999.0)
        encoding = {}
        for var in ds_out.data_vars:
            encoding[var] = {'zlib': True, 'complevel': 4, '_FillValue': -9999.0}
        ds_out.to_netcdf(outfile_ens, encoding=encoding)

    t2 = time.time()
    print(f'Complete prob. estimation (ensemble generation) of {target_vars} with linked vars {linkvar} for member range {member_range_2}.')
    print(f'Time cost (s): {t2-t1}')

def generate_prob_estimates(config):
    t1 = time.time()

    ensemble_start = config['ensemble_start']
    ensemble_end   = config['ensemble_end']
    num_processes  = config['num_processes']

    print('#' * 50)
    print('Probabilistic estimation')
    print('#' * 50)
    print('Target variables:', config['target_vars'])
    print('Number of processes:', num_processes)

    if ensemble_start<0 or ensemble_end<0:
        print('Error: ensemble_start or ensemble_start is negative.')
        return

    spcorr_structure(config)

    ensemble_number = ensemble_end - ensemble_start + 1
    items = [(config, [e, e+1]) for e in range(ensemble_number)]
    with Pool(num_processes) as pool:
        pool.starmap(generate_prob_estimates_serial, items)

    t2 = time.time()
    print('Time cost (s):', t2 - t1)
    print('Probabilistic estimation (ensemble generation) completed successfully!\n\n')

