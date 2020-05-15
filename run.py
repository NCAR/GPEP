import numpy as np
import os

def correct_merge(stndata, readata_raw, readata_stn, reacorr_stn, reamerge_stn, reamerge_weight_stn, neargrid_loc,
         neargrid_dist, merge_choice, mask, var, hwsize, corrmode, anombound):

    nrows, ncols, nearnum = np.shape(neargrid_loc)
    reanum, nstn, nday = np.shape(readata_stn)

    # correct raw gridded reanalysis data using all stations
    readata_corr = np.nan * np.zeros([reanum, nrows, ncols, nday], dtype=np.float32)
    for rr in range(reanum):
        # calculate correction ratio at all station point
        anom_ori = calculate_anomaly(readata_stn[rr, :, :], stndata[:, :],
                                     hwsize, corrmode, upbound=anombound[1], lowbound=anombound[0])
        anom_ext = extrapolation(anom_ori, neargrid_loc, neargrid_dist)
        readata_corr[rr, :, :, :] = error_correction(readata_raw, anom_ext, mode=corrmode)

    # merge reanalysis data
    merge_data = np.nan * np.zeros([nrows, ncols, nday], dtype=np.float32)
    weight_grid = extrapolation(reamerge_weight_stn, neargrid_loc, neargrid_dist)
    for i in range(nday):
        datain = np.zeros([nrows, ncols, reanum])
        for rr in range(reanum):
            datain[:, :, rr] = readata_corr[rr, :, :, i]
        merge_data[:, :, i] = weightmerge(datain, weight_grid)

    # calculate the error of merged data (this is actually independent with merged data estimation)
    error1 = extrapolation(reamerge_stn - stndata, neargrid_loc, neargrid_dist)
    error2 = np.nan * np.zeros([reanum, nrows, ncols, nday])
    for rr in range(reanum):
        error2[rr, :, :, :] = extrapolation(reacorr_stn[rr, :, :] - stndata, neargrid_loc, neargrid_dist)

    merge_error = np.nan * np.zeros([nrows, ncols, nday], dtype=np.float32)
    for r in range(nrows):
        for c in range(ncols):
            if not np.isnan(mask[r, c]):
                chi = merge_choice[r, c]
                if chi > 0:
                    merge_error[r, c, nday] = error2[chi - 1, r, c, :]
                else:
                    merge_error[r, c, nday] = error1[r, c, :]

    return readata_corr, merge_data, merge_error