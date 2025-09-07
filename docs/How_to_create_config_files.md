
# How to create configuration files for GPEP

To ensure that your configuration files are correct and efficient, it is recommended that you base them on the template files in the `./config_templates` folder. The format for these files is [TOML](https://toml.io/en/), which is easy to read and simple to use.

When running a specific case, you will need two configuration files: the control file and the setting file. These files can be found in `./config_templates` under the names `testcase.config.dynamic.toml` and `model.settings.toml`, respectively. The file name of `model.settings.toml` is defined in `testcase.config.dynamic.toml` to ensure that the model knows what settings to use. The setting file only needs minimal or no changes unless major data or method choice differences exist between a specific case and the test case. Additionally, there are other example configuration files, such as `testcase.config.RF.toml`, which provide alternate settings.  

Besides, to create a test case from scratch, we recommend users follow the test case to prepare all necessary files. The test cases of GPEP are open access on [Zenodo](https://zenodo.org/record/8222852). They can be obtained using the `./tools/get_testcase.py` script, which is used in `./docs/GPEP_demo.ipynb`. See `./tools/README.md` for more descriptions. The example files shown in this document are from the test case `../GPEP_test_cases/cali2017` which is the default folder in `./tools/get_testcase.py`. The `../` relative path assumes you are in the main directory of GPEP.

For a more detailed explanation, please refer to the following section.

# Configuration file
This section explains the parameters of the configuration file `testcase.config.dynamic.toml`. Please note that the file name can be changed by users.

## case_name

This parameter is used to name GPEP outputs.

## num_processes

This parameter determines the number of workers or CPU cores used for parallelization. If generating ensemble outputs, different members will run in parallel. For example, if you want to generate 10 members but assign 16 cores to the job, the remaining 6 cores will be idle. However, all 16 cores will be used in the spatial regression part, which could be more time-consuming than ensemble generation.

## modelsettings_file

This parameter specifies the GPEP model settings file, e.g., `model.settings.toml` in the `../GPEP_test_cases/cali2017`. This file does not need to be changed for a specific case using default settings. Still, it contains crucial and advanced parameters (e.g., machine learning methods).

## input_stn_list
*Unnecessary if input_stn_all is provided*  

This parameter specifies the CSV file containing the IDs and attributes of stations used by GPEP to generate gridded estimates. The CSV file should have the `stnid`, `lat`, and `lon` columns, and optionally, other columns with auxiliary variables used in the regression. If you would like to use variables like elevation and slope, make sure you include those variables in this file and define the variable names in the configuration file. An example file is `../GPEP_test_cases/cali2017/inputs/CALI.screened_stn_list_slope.v3.csv`.

Below is a simplified example table of a station list.  
| stnid       | lat     | lon        | elev  | slp_n    | slp_e  |
| ----------- | ------- | ---------- | ----- | -------- | ------ |
| US1CAAL0003 | 37.7169 | \-122.0585 | 163.4 | 1.953    | 4.542  |
| US1CADN0001 | 41.824  | \-124.1482 | 19.8  | \-2.948  | 15.123 |
| US1CAFR0003 | 36.9897 | \-119.3971 | 509.9 | \-18.502 | 27.648 |
| US1CAHM0001 | 40.8809 | \-124.0692 | 128   | 1.367    | 13.23  |
| US1CAHM0002 | 40.8248 | \-124.0382 | 23.5  | 1.73     | 14.267 |
| US1CAHM0003 | 41.2084 | \-124.1027 | 172.8 | \-0.565  | 12.834 |
| US1CAHM0004 | 40.9231 | \-124.084  | 99.1  | 1.367    | 13.23  |

## input_stn_path

*Unnecessary if input_stn_all is provided*  

This parameter specifies the path of station inputs. The file names are `stnid.nc`, where `stnid` comes from `input_stn_list`. The input variables (e.g., `prcp`, `tmin`, and `tax`) are read from these stations. Other auxiliary variables (e.g., `fill_flag_prcp`) will not be used by GPEP. Note that it is okay if some input variables are missing in some station files because sometimes stations may only measure a few variables. An example folder path is `../GPEP_test_cases/cali2017/inputs/stndata`, and below is a simplified example of a station file named `USW00094299.nc` .
```
netcdf USW00094299 {
dimensions:
	time = UNLIMITED ; // (28 currently)
variables:
	int fill_flag_prcp(time) ;
		fill_flag_prcp:long_name = "Fill flag" ;
	int fill_flag_tmax(time) ;
		fill_flag_tmax:long_name = "Fill flag" ;
	int fill_flag_tmin(time) ;
		fill_flag_tmin:long_name = "Fill flag" ;
	float prcp(time) ;
		prcp:long_name = "Total daily precipitation" ;
	double time(time) ;
		time:calendar = "standard" ;
		time:units = "seconds since 1970-01-01 00:00:00" ;
	float tmax(time) ;
		tmax:long_name = "Maximum daily air temperature" ;
	float tmin(time) ;
		tmin:long_name = "Minimum daily air temperature" ;
}
```

## input_stn_all

The `input_stn_all` parameter is used to provide station input data to GPEP. It is not necessary to use this parameter if `input_stn_list` and `input_stn_path` are already provided. However, if all three parameters are provided, `input_stn_all` has a higher priority.

This parameter is particularly useful when all stations and variables are stored in a single netCDF file, which can be more convenient in some cases. For example, the file `../GPEP_test_cases/cali2017/outputs/stn_info/all_station_data.nc` is generated by GPEP using `input_stn_list` and `input_stn_path`, and can be used as a template for creating similar files for other domains. To use `input_stn_all`, simply provide the path to the netCDF file containing all station and variable data.
```
netcdf all_stn {
dimensions:
	stn = 1219 ;
	time = 28 ;
	string11 = 11 ;
variables:
	int64 stn(stn) ;
	int64 time(time) ;
	char stnid(stn, string11) ;
	double lat(stn) ;
	double lon(stn) ;
	double elev(stn) ;
	double slp_n(stn) ;
	double slp_e(stn) ;
	float prcp(stn, time) ;
	float tmin(stn, time) ;
	float tmax(stn, time) ;
	float tmean(stn, time) ;
	float trange(stn, time) ;
}
```

## infile_grid_domain
This netCDF file defines the target domain. An example can be found below. Note that its dimensions are defined as `x` and `y` instead of `lat` and `lon`. The grids, such as latitude and longitude, can be irregular.

```
netcdf CALI.gridinfo.0625_v3 {
dimensions:
	param = 1 ;
	y = 144 ;
	x = 112 ;
variables:
	double dx(param) ;
		dx:long_name = "Grid spacing in x-direction" ;
		dx:units = "degrees" ;
	double dy(param) ;
		dy:long_name = "Grid spacing in y-direction" ;
		dy:units = "degrees" ;
	double elev(y, x) ;
		elev:name = "Elevation" ;
		elev:long_name = "Elevation of grid" ;
		elev:units = "m" ;
		elev:_FillValue = -999. ;
	double gradient_n_s(y, x) ;
		gradient_n_s:name = "N_S_Slope" ;
		gradient_n_s:long_name = "Smoothed Topographic Gradient (Slope) from North to South" ;
		gradient_n_s:units = "-" ;
		gradient_n_s:_FillValue = -999. ;
	double gradient_w_e(y, x) ;
		gradient_w_e:name = "W_E_Slope" ;
		gradient_w_e:long_name = "Smoothed Topographic Gradient (Slope) from West to East" ;
		gradient_w_e:units = "-" ;
		gradient_w_e:_FillValue = -999. ;
	double latitude(y, x) ;
		latitude:long_name = "Latitude" ;
		latitude:units = "degrees_north" ;
	double longitude(y, x) ;
		longitude:long_name = "Longitude" ;
		longitude:units = "degrees_east" ;
	int mask(y, x) ;
		mask:long_name = "Basin mask for running ensemble" ;
		mask:units = "-" ;
	double startx(param) ;
		startx:long_name = "Lower left longitude" ;
	double starty(param) ;
		starty:units = "degrees" ;
		starty:long_name = "Lower left latitude" ;
```


## **outpath_parent**

The `outpath_parent` parameter specifies the folder where all GPEP outputs will be saved.

## **dynamic_predictor_filelist**

_Optional._ The `dynamic_predictor_filelist` parameter is a list of input files that contain dynamic predictors used to support gridding estimation. If the file cannot be found (e.g., "NA" file name) or is not provided, dynamic predictors will not be used in the estimation process. The dynamic input netcdf files must have structured/regular latitude and longitude dimensions. The dynamic inputs will be interpolated to station points and regular/irregular grids of the target domain.

An example of the `dynamic_predictor_filelist` parameter is provided in `../GPEP_test_cases/cali2017/inputs/gridded_datafile_list.txt`, which is shown below:

```
../GPEP_test_cases/cali2017/inputs/griddata/inputfile_0.nc
../GPEP_test_cases/cali2017/inputs/griddata/inputfile_1.nc
```

Part of the information of `inputfile_0.nc` is as below:
```
netcdf inputfile_0 {
dimensions:
	time = 5 ;
	lat = 144 ;
	lon = 112 ;
variables:
	double cube_root_prec_rate(time, lat, lon) ;
		cube_root_prec_rate:long_name = "Precipitation rate" ;
	double tmp_2m(time, lat, lon) ;
		tmp_2m:long_name = "Temperature" ;
	int64 time(time) ;
	double lon(lon) ;
	double lat(lat) ;
}
```


## date_start

Start date of the target output period in YYYY-MM-DD format.

## date_end

End date of the target output period in YYYY-MM-DD format. Note that `date_start` and `date_end` should be a subset of the period of raw input stations (`input_stn_list`). The raw station records can contain a longer period of data.

## input_vars

A list of input variables to be read from raw station files, such as `['prcp', 'tmin', 'tmax']`.

## target_vars

A list of target output variables. This does not necessarily have to be the same as the `input_vars`.

## mapping_InOut_var

_Optional._ A list that maps input variables to target variables. If the input variable and target variable have the same name, there is no need to add a mapping relation. Examples include:

-   Converting tmin and tmax to tmean and trange: `mapping_InOut_var = ['tmean=(tmin+tmax)/2', 'trange=abs(tmax-tmin)']`
-   Converting precipitation unit from inch to cm: `mapping_InOut_var = ['prcp=prcp * 2.54']`

## target_vars_WithProbability

_Optional._ A list of variables (e.g., precipitation, snowfall, hail) requiring occurrence estimation with a probability between 0-1. This list should be a subset of `target_vars`, such as `['prcp']`.  

## probability_thresholds

_Optional._ Thresholds used to define the occurrence of variables from target_vars_WithProbability. Default values are zero for all variables if probability_thresholds is not provided or empty. For example, if `target_vars_WithProbability=['prcp']`, you can use `probability_thresholds=[0]`. Higher thresholds can be used to focus on large events.

## minRange_vars

_Optional._ A list corresponding to `target_vars` that defines the minimum values of the corresponding variable in spatial regression. For variables with minimum values such as 0 for precipitation, this parameter is useful because regression models could get negative values. If not provided, no constraint will be applied. Setting the value to `-inf` can also turn off the constraint. For example, `minRange_vars = [0, -100, 1]`.

## maxRange_vars

_Optional._ Similar to `minRange_vars` but for maximum value constraint. Setting the value to `inf` or `+inf` can also turn off the constraint. For example, `maxRange_vars = [2000, 60, 60]`.

## transform_vars

_Optional._ Some variables (e.g., precipitation) may need to be transformed to reduce the skewness of the probability distribution. `boxcox` transformation is supported with the parameter defined in the `model.settings.toml` file. An empty value means no transformation. If not provided, no transformation. For example, `transform_vars = ['boxcox', '', '']`. For probabilistic estimation of skewed variables such as precipitation, transformation is recommended. In contrast, if only deterministic estimation is needed, transformation is not necessary.

## dynamic_predictor_name

_Optional._ A list corresponding to `target_vars` that defines what dynamic predictors will be used in gridding estimation for each target variable. Invalid variable names or an empty list means that dynamic predictors won't be used. For example, `dynamic_predictor_name = ['cube_root_prec_rate, tmp_2m', 'tmp_2m', '']` where `prcp` will use `cube_root_prec_rate` and `tmp_2m` as dynamic predictors, `tmean` will use `tmp_2m`, and `trange` will not use dynamic predictors.

```
dynamic_predictor_name = [['cube_root_prec_rate', 'tmp_2m'],  
                          ['tmp_2m'],  
                          [],  
                          ]
```

## predictor_name_static_stn
This list contains names of static predictors of input stations that will be used to support gridding estimation. These variables must be contained in `input_stn_list` or `input_stn_all`.  
Example: `predictor_name_static_stn = ['lat', 'lon', 'elev', 'slp_n', 'slp_e']`

## predictor_name_static_grid

This list corresponds to `predictor_name_static_stn` but are static predictors of grids contained in `infile_grid_domain`. For example, `predictor_name_static_grid = ['latitude', 'longitude', 'elev', 'gradient_n_s', 'gradient_w_e']`


## dynamic_predictor_operation

_Optional._ Dynamic predictors may need some processing, such as interpolation from dynamic grids to target grids and station points, or transformation of the raw dynamic predictor.

The keyword `interp` accepts the following methods, enabled by `xarray.interp`: {"linear", "nearest", "zero", "slinear", "quadratic", "cubic", "polynomial"}. The default interpolation method is `nearest`.

The keyword `transform` accepts methods defined in the transformation section.

Example: `cube_root_prec_rate:interp=linear:transform=boxcox` or `dynamic_predictor_operation = []`. If the `transform` keyword is missing, no transformation will be performed.

## ensemble_flag

Values: `true` or `false`.

If `true`, GPEP will generate ensemble members of probabilistic estimates. If `false`, GPEP will not generate probabilistic estimates and will stop after generating spatial gridded estimates which are deterministic, i.e., GPEP will act as an interpolation tool. In this case, the below ensemble settings will not play a part in GPEP.



## ensemble_start

The start member of ensemble outputs. It starts from 1, not 0.

## ensemble_end

The end member of ensemble outputs. It starts from 1, not 0.


## linkvar

For GMET probabilistic estimation, some variables (e.g., precipitation) depend on other variables (e.g., temperature range) when generating spatiotemporally correlated random fields.

For example, for `linkvar = [ ['prcp', 'trange'], ]`, `prcp` depends on `trange`. Note, `trange` does not depend on `prcp` in this setting. To turn it off, use `linkvar = []`.

## clen

_Optional._ This is a list corresponding to `target_vars`. It is the spatial correlation length (km). If not provided or the provided values are negative, GPEP will estimate `clen` from raw input station data. Note that GPEP will use all input station data (i.e., `input_stn_list`) instead of the target period between `date_start` and `date_end` which can be short.

## lag1_auto_cc

Similar to `clen`, but for the lag-1 autocorrelation of a variable. If not provided or the provided values are smaller than -1, GPEP will estimate it from raw station inputs.

## cross_cc

This is a list corresponding to `linkvar`, not `target_vars`. It is the cross-correlation between the time series of two variables. The parameters related to cross-correlation calculation are in `model.settings.toml`.

## auto_corr_method

`direct`: Using raw time series to calculate autocorrelation. `anomaly`: Uses an `n`-timestep (e.g., 31-day) moving average (window) to remove temporal cycle. For example, `rolling_window=31`

## target_vars_max_constrain
*optional*. For certain variables, such as precipitation, the maximum values (with scale factors) from nearby weather stations can be utilized to constrain probabilistic estimates. It is essential to exercise caution when using this method because it may not always be appropriate.

If the `target_vars_max_constrain` list is empty or not provided, the maximum value constraint will not be applied.

Example: `target_vars_max_constrain = ['prcp']`  

# Model setting file
Example: `model.settings.toml`. The file name can be changed and should be defined as `modelsettings_file` in the configuration file.  

## master_seed

This parameter controls the random seed used for probabilistic processes such as machine learning and probabilistic estimation. If the value is negative, the generation is truly random without reproducibility.

## gridcore_continuous
Gridding method for continuous spatial regression.  
Sklearn module is used to support most functions: https://scikit-learn.org/stable/supervised_learning.html
- Locally weighted regression.
(1) The original GMET method is `LWR:Linear`.
(2) Sklearn-based methods support simple usage with `model.fit()` and `model.predict` or `model.predict_prob`, in the format of `LWR:linear_model.METHOD`. Examples of METHOD are `LinearRegression`, `LogisticRegression`, `Ridge`, `BayesianRidge`, `ARDRegression`, `Lasso`, `ElasticNet`, `Lars`, etc.  

- Global regression using machine learnig methods:
Machine learning methods are supported by sklearn. Parametrs of methods supported by sklearn can be defined at the bottom of this configuration file (optional). Examples: `ensemble.RandomForestRegressor`, `ensemble.RandomForestClassifier`, `neural_network.MLPRegressor`, `neural_network.MLPClassifier`, `ensemble.GradientBoostingClassifier`, and `ensemble.GradientBoostingRegressor`.

The parameters of sklearn methods can be defined in the [sklearn] section.  

Example: `gridcore_continuous = 'LWR:Linear'`

## gridcore_classification
Same with `gridcore_continuous` but for estimating the probability of event.  
Example: `gridcore_continuous = 'LWR:Logistic'`


## n_splits

This parameter specifies the number of batches in which all samples are divided for cross-validation in machine learning-based global regression methods. For local regression, leave-one-out cross-validation is implemented.

## nearstn_min and nearstn_max

These parameters define the minimum and maximum numbers of nearby stations to be used for regression of a target point.

## try_radius

When searching for nearby stations, this parameter defines the initial search radius (km). If not enough stations are found, the radius is expanded.

## initial_distance
Initial Search Distance in km (expanded if need be).  
Example: `initial_distance = 100`

## weight_formula
Weight calculation formula. Only two variables/parameters are allowed in the formula.   
- dist: distance between points (km). For example,   
- maxdist (optional): max(initial_distance, max(dist)+1), which is a parameter used in weight calculation   

The default formula is `weight_formula = '(1 - (dist / maxdist) ** 3) ** 3'` following Python syntax. 3 is the exponential factor and is the default parameter. Users can change the formula or the exponential factor.  

## [netcdf_info]

### stn_lat_name and stn_lon_name

These parameters define the names of the latitude and longitude dimensions in the input station files.

### grid_lat_name and grid_lon_name

These parameters define the names of the latitude and longitude dimensions in the 2D array of the target domain.

### grid_mask_name

This parameter defines the name of the 2D mask array of the target domain. A value of 1 indicates that the grid will be generated, while a value of 0 indicates that it will not.

### dynamic_grid_lat_name and dynamic_grid_lon_name

These parameters define the names of the latitude and longitude dimensions in the dynamic predictor files.


## [transform]
Transformation section. This will be read as a dictionary in Python. Currently,`boxcox` and `ecdf` transformation is supported. Please note that `ecdf` is still experimental, and requires a longer station data record to be provided to build the empirical cdf (ecdf).

### [transform.boxcox]
#### exponent
The exponential factor used in box-cox transformation. Example: `exponent = 4`.

### [transform.ecdf]
#### pooled
Boolean to describe if a pooled (True) or station ecdf approach is taken. Example: `pooled=True`.
#### min_z_value
This is a value to be assigned to all nan values, following the ecdf transform. It represents a lower standard score (z), in the case of precipitation, it corresponds to a near zero amount. Example: `min_z_value=-4`.
#### interp_method
Selection of the method to translate cumulative probabilities back to measurements values in the ecdf back transform. Current options are linear interpolation with extrapolation, and a gamma distribution fit.
Example: `interp_method=interp1d`, `interp_method=gamma`.

#### min_est_value
Very small values from back transformation below this threshold value are converted to nan, and then zero. Minimum estimated value is used a a cut off for mariginal values that are created in the regression process.
Example: `min_est_value=0.01`.


## [sklearn]
*Optional*
scikit-learn section defining the parameters of sklearn methods used by `gridcore_continuous` and `gridcore_classification`. If no parameters are provided or if the section does not even exist, default parameters will be used.    
Note: Use the pure name of sklearn method as the sub-dict name. For example, if `gridcore_continuous='ensemble.RandomForestRegressor'`, define the parameters of `ensemble.RandomForestRegressor` under the section `[sklearn.RandomForestRegressor]` (see the example). Don't inlcude `ensemble.` in the dict name. Similarly, if `gridcore_continuous='LWR:linear_model.Ridge'`, define the parameters of `Ridge` under the section `sklearn.Ridge`.

Example:
### [sklearn.RandomForestRegressor]
#### n_estimators   = 500
#### n_jobs = 5

## [flags]
Boolean flags. They don't affect model performance. For simple run, no need to change flag values.  

### output_random_fields

If set to `true`, random fields will be output in the ensemble outputs. If set to `false` (default), random fields will not be output to reduce file sizes.

### append_date_to_output_filename

If set to `true` (default), output file names will contain the date range. If set to `false`, output file names will not contain date range.

### print_config
If set to `true`, all configurations will be printed at the beginning of model printouts. If set to `false` (default), nothing will be printed for simplicity.  

### overwrite flags
**overwrite_stninfo**: Files related integrated station inputs and nearby information
**overwrite_station_cc**: For the station correlation file  
**overwrite_weight**: For station weight file  
**overwrite_cv_reg**: For station-based cross validation file  
**overwrite_grid_reg**: For gridded regression file  
**overwrite_ens**: For ensemble output files  
**overwrite_spcorr**: For spatial correlation structure files  
By default, these flags are all `false` to utilize existing outputs to speed up production. Users may change them to `true` to generate brand new files for all or part of GPEP outputs to overwrite existing files. A safer way is to just change `outpath_parent` to create files in a new folder.  
