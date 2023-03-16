import read_config
import data_processing
import near_stn_search
import weight_calculation
import regression
import probabilistic_auxiliary
import data_correlation
import probabilistic_estimation

import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':

    config_file = "../test_cases/testcase.config_short.toml"
    # config_file = sys.argv[1]

    ########################################################################################################################
    # load configuration file

    config = read_config.read_config(config_file)

    ########################################################################################################################
    # assemble individual stations and station attributes (e.g., lat, lon) to one netcdf file
    config = data_processing.assemble_fortran_GMET_stns_to_one_file(config)

    ########################################################################################################################
    # get near station info for each station/grid
    config = near_stn_search.get_near_station_info(config)

    ########################################################################################################################
    # calculate weights based on near station info. this step is independent to enable flexible weight test if needed.
    config = weight_calculation.calculate_weight_using_nearstn_info(config)

    ########################################################################################################################
    # regression:
    # cross validation station regression (i.e., at station points)
    config = regression.main_regression(config, 'loo')

    # grid regression
    config = regression.main_regression(config, 'grid')

    ########################################################################################################################
    # probabilistic / ensemble estimation
    if config['ensemble_flag'] == False:
        print('ensemble_flag if false in the configuration file. No need to generate ensemble estimation.')
    else:
        ########################################################################################################################
        # get uncertainty estimation based on difference between cross validation station regression estimates and station observations
        # interpolation from points to grids
        config = probabilistic_auxiliary.extrapolate_auxiliary_info((config))

        ########################################################################################################################
        # probabilistic estimation

        # 1. get space and time correlations
        config = data_correlation.station_space_time_correlation(config)

        # 2. probabilistic estimation
        config = probabilistic_estimation.generate_probabilistic_estimates(config)

    print('Successfully finish the program!')

