import read_config
import data_processing
import near_stn_search
import weight_calculation
import regression
import probabilistic_auxiliary
import data_correlation
import probabilistic_estimation

import sys, time
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
    t1 = time.time()

    config_file = sys.argv[1]

    ########################################################################################################################
    # load configuration file

    config = read_config.read_config(config_file)

    ########################################################################################################################
    # assemble individual stations and station attributes (e.g., lat, lon) to one netcdf file
    config = data_processing.merge_stndata_into_single_file(config)

    ########################################################################################################################
    # get near station info for each station/grid
    config = near_stn_search.get_near_station_info(config)

    ########################################################################################################################
    # calculate weights based on near station info. this step is independent to enable flexible weight test if needed.
    config = weight_calculation.calculate_weight_using_nearstn_info(config)

    ########################################################################################################################
    # perform regression
    # (1) estimate predictive uncertainty using cross-validated (i.e., leave one out, LOO) regression at station points 
    config = regression.main_regression(config, 'cval')

    # (2) estimate regression coefficients at all grid points
    config = regression.main_regression(config, 'grid')

    ########################################################################################################################
    # probabilistic / ensemble estimation
    if config['ensemble_flag'] == False:
        print('ensemble_flag is false in the configuration file -- ensemble generation will be skipped.')
        
    else:
        ########################################################################################################################
        # estimate gridcell uncertainty based on interpolating estimation error from station locations (the LOO regression)
        config = probabilistic_auxiliary.extrapolate_auxiliary_info((config))

        ########################################################################################################################
        # probabilistic estimation (ensemble generation)

        # 1. get space and time correlations
        config = data_correlation.station_space_time_correlation(config)

        # 2. probabilistic estimation
        config = probabilistic_estimation.generate_prob_estimates(config)

    t2 = time.time()
    print('Total time cost (s):', t2-t1)
    print('Successfully finished PyGMET run!')

