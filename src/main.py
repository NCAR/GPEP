import glob

import toml, os, time, sys
import pandas as pd
import numpy as np
import xarray as xr

import data_processing
import near_stn_info

config_file = "example.PyGMET.config.toml"

########################################################################################################################
# load configuration file

config = toml.load(config_file)
print(config)


########################################################################################################################
# assemble individual stations and station attributes (e.g., lat, lon) to one netcdf file
config = data_processing.assemble_fortran_GMET_stns_to_one_file(config)

########################################################################################################################
# assemble individual stations and station attributes (e.g., lat, lon) to one netcdf file
config = near_stn_info.get_near_station_info_and_weight(config)


########################################################################################################################
# 