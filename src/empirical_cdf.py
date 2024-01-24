# calculation and transformation of data via empirical CDFs

import os, time, sys
import pandas as pd
import numpy as np
import xarray as xr

from scipy.stats import norm
from scipy.interpolate import interp1d

def calculate_monthly_cdfs(ds,var_name):
    """
    Calculate the empirical cumulative distribution functions (CDFs) for each month and each station.

    This function computes the empirical CDF for the precipitation data of each station for each month,
    capturing the distribution of daily precipitation values over the years.

    :param all_stn_path: Path to the netCDF file containing the precipitation data for all stations.

    :return: A nested dictionary where the top-level keys are station identifiers, and each value is another
             dictionary with months as keys (1 through 12) and the corresponding CDF data as values.
    """

    # Read in the data
    print(f'Reading data from {var_name}')
    df = ds[var_name].to_dataframe()

    all_stations_cdfs = {}
    for station in df.index.get_level_values(0).unique():
        monthly_cdfs = {}
        for month in range(1, 13):
            # Extract data for the specific station
            df_station = df.loc[station]

            # Extract data for the specific month and station, dropping missing values
            month_data = df_station[df_station.index.month == month]

            #Remove any zero or negative values
            month_data = month_data[month_data > 0]
            month_data_w_nan = month_data.values
            month_data_values = month_data_w_nan[~np.isnan(month_data_w_nan)]

            # Sort the data and calculate the empirical CDF values
            sorted_data = np.sort(month_data_values)
            cdf_values = np.arange(1, len(sorted_data) + 1) / len(sorted_data)

            # Store the CDF for this month and station
            monthly_cdfs[month] = pd.DataFrame({'Value': sorted_data, 'CDF': cdf_values})

        # Store the CDFs for all months for this station
        all_stations_cdfs[station] = monthly_cdfs

    return all_stations_cdfs

def normal_quantile_transform(df, all_stations_cdfs):
    """
    Apply a normal quantile transform to the precipitation data of each station for each month.

    This function matches the empirical cumulative probabilities of precipitation values
    with the cumulative probabilities from a standard normal distribution.

    :param df: DataFrame containing precipitation data, indexed by date, with columns as stations.
    :param all_stations_cdfs: Dictionary containing the empirical CDFs for each station and month.
    :return: DataFrame containing the Z-scores after transformation.
    """
    transformed_data = pd.DataFrame(index=df.index)

    for station in df.columns:
        if sum(df[station]) == 0:
            print(f'No precipitation data for station {station}, skipping station in normal_quantile_transform()')
            continue
        for month in range(1, 13):
            month_data = df[station][df.index.month == month]

            #Remove any zero or negative values
            month_data = month_data[month_data > 0]

            empirical_cdf = all_stations_cdfs[station][month]

            if empirical_cdf is not None and not empirical_cdf.empty:
                cdf_interp = interp1d(empirical_cdf['Value'], empirical_cdf['CDF'],bounds_error=False, fill_value=(-2.55,2.55))
                cum_probs = np.clip(cdf_interp(month_data.dropna()), 0, 0.999)
                z_scores = norm.ppf(cum_probs)
                transformed_data.loc[month_data.index, station] = z_scores
            else:
                transformed_data.loc[month_data.index, station] = np.nan
                #print(f'Missing data for empirical_cdf for station {station}, month {month}')

    return transformed_data

def inverse_normal_quantile_transform(df, all_stations_cdfs):
    """
    Reverse the normal quantile transform applied to the precipitation data.

    This function maps the Z-scores back to the original precipitation values using
    the empirical cumulative distribution functions (CDFs) for each station and month.

    :param transformed_data: DataFrame containing the Z-scores, indexed by date, with columns as stations.
    :param all_stations_cdfs: Dictionary containing the empirical CDFs for each station and month.
    :return: DataFrame containing the original precipitation values after reverse transformation.
    """
    back_transformed_data = pd.DataFrame(index=df.index)

    for station in df.columns:
        for month in range(1, 13):

            z_scores = df[station][df.index.month == month]

            empirical_cdf = all_stations_cdfs[station][month]

            if empirical_cdf is not None and not empirical_cdf.empty:

                cum_probs = norm.cdf(z_scores)

                value_interp = interp1d(empirical_cdf['CDF'], empirical_cdf['Value'], bounds_error=False, fill_value=(0,100))

                original_values = value_interp(cum_probs)

                #Fill nan values with zero
                original_values = np.nan_to_num(original_values)

                back_transformed_data.loc[z_scores.index, station] = original_values
            else:
                # Handle empty or missing data appropriately (e.g., set to NaN)
                back_transformed_data.loc[z_scores.index, station] = np.nan

    back_transformed_data = back_transformed_data.T
    
    return back_transformed_data