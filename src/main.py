import glob

import toml, os, time, sys
import pandas as pd
import numpy as np
import xarray as xr

import data_processing
import near_stn_search
import weight_calculation
import regression
import spatial_extrapolation
import data_correlation

config_file = "example.PyGMET.config.toml"

########################################################################################################################
# load configuration file

config = toml.load(config_file)
print(config)


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
# leave-one-out station regression (i.e., at station points)
config = regression.main_regression(config, 'loo')

# grid regression
config = regression.main_regression(config, 'grid')

########################################################################################################################
# get uncertainty estimation based on difference between leave-one-out station regression estimates and station observations
# interpolation from points to grids
config = spatial_extrapolation.station_error_extrapolation((config))


########################################################################################################################
# probabilistic estimation
config = data_correlation.station_space_time_correlation(config)




