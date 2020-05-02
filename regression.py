import numpy as np
import auxiliary as au
import sys

def least_squares(x, y, tx):
    # In fortran version, ludcmp and lubksb are used to calcualte matrix inversion
    # call ludcmp(a, indx, d)
    # call lubksb(a, indx, b)

    # In Python version, numpy is used to calculate matrix inversion
    b = np.matmul(tx, y)
    a = np.matmul(tx, x)

    deta = np.linalg.det(a)  # Compute the determinant of an array
    if deta == 0:
        print('Singular matrix')
        b[:] = 0
    else:
        ainv = np.linalg.inv(a)
        b = np.matmul(ainv, b)

    return b


def logistic_regression(x, tx, yp):
    # nstn: station number
    # nvars: station attributes (1, lat/lon/...), 1 is for regression
    nstn, nvars = np.shape(x)

    b = np.zeros(nvars)  # regression coefficients (beta)
    p = np.zeros(nstn)  # estimated pop (probability of precipitation)
    f = 0  # flag: continue or stop loops
    it = 0  # iteration times

    while f != 1:
        # check for divergence
        xb = -np.matmul(x, b)
        if np.any(xb > 50):
            f = 1
        else:
            p = 1 / (1 + np.exp(xb))

        # check for divergence
        if np.any(p > 0.9999):
            # logistic regression diverging
            f = 1
        else:
            v = np.zeros([nstn, nstn])  # diagonal variance matrix
            for i in range(nstn):
                v[i, i] = p[i] * (1 - p[i])
            xv = np.matmul(v, x)
            yn = yp - p  # difference between station occurrence (0/1) and estimated pop: Bnew -Bold in Martyn and Slater 2006
            bn = least_squares(xv, yn, tx)

            # check: converging
            if np.any(np.abs(bn) > 1e-4):
                f = 0
            else:
                f = 1

            # check: iteration times
            if it > 20:
                f = 1
            else:
                f = 0

            b = b + bn  # update coefficients

        it = it + 1

    return b


def station_error(prcp_stn, tmean_stn, trange_stn, stninfo, near_stn_prcpLoc, near_stn_prcpWeight, near_stn_tempLoc,
                  near_stn_tempWeight, trans_exp, trans_mode, nearstn_min):
    nstn, ntimes = np.shape(prcp_stn)
    pcp_err_stn = -999 * np.ones([nstn, ntimes])
    tmean_err_stn = -999 * np.ones([nstn, ntimes])
    trange_err_stn = -999 * np.ones([nstn, ntimes])

    for t in range(ntimes):
        print('Current time:', t, 'Total times:', ntimes)
        # assign vectors of station alues for prcp_stn, temp, for current time step
        # transform prcp_stn to approximate normal distribution
        y_prcp = au.transform(prcp_stn[:, t], trans_exp, trans_mode)
        y_tmean = tmean_stn[:, t]
        y_trange = trange_stn[:, t]

        for gg in range(nstn):
            if prcp_stn[gg, t] > -1:
                # reduced matrices for precip
                nstn_prcp = int(np.sum(near_stn_prcpLoc[gg, :] > -1))
                if nstn_prcp >= nearstn_min:
                    w_pcp_red = np.zeros([nstn_prcp, nstn_prcp])
                    for i in range(nstn_prcp):
                        w_pcp_red[i, i] = near_stn_prcpWeight[gg, i]  # eye matrix: stn weight in one-one line
                    w_pcp_1d = near_stn_prcpWeight[gg, 0:nstn_prcp]  # stn weight
                    w_pcp_1d_loc = near_stn_prcpLoc[gg, 0:nstn_prcp]  # stn ID number/location
                    y_prcp_red = y_prcp[w_pcp_1d_loc]  # transformed prcp_stn
                    x_red = stninfo[w_pcp_1d_loc, :]  # station lat/lon/ele/slope_ns/slope_we

                    yp_red = np.zeros(nstn_prcp)  # pop: 0/1
                    yp_red[prcp_stn[w_pcp_1d_loc, t] > 0] = 1
                    ndata = np.sum(yp_red == 1)  # number of prcp_stn>0
                else:
                    # there are not enough nearby stations
                    x_red = 0  # not really necessary. just to stop warming from Pycharm
                    w_pcp_red = 0  # not really necessary. just to stop warming from Pycharm
                    yp_red = 0  # not really necessary. just to stop warming from Pycharm
                    y_prcp_red = 0  # not really necessary. just to stop warming from Pycharm
                    ndata = 0

                # prcp processing
                if ndata == 0:
                    # nearby stations do not have positive prcp data
                    pcp_err_stn[gg, t] = 0
                else:
                    # tmp needs to be matmul(TX, X) where TX = TWX_red and X = X_red
                    mat_test = np.matmul(np.transpose(x_red), w_pcp_red)
                    tmp = np.matmul(mat_test, x_red)
                    vv = np.max(np.abs(tmp), axis=1)

                    # decide if slope is to be used in regression
                    if np.any(vv) == 0 or (abs(stninfo[gg, 4]) < 3.6 or abs(stninfo[gg, 5]) < 3.6):
                        slope_flag_pcp = 0
                        x_red_use = x_red[:, 0:4]  # do not use slope in regression
                        stninfo_use = stninfo[gg, 0:4]
                    else:
                        slope_flag_pcp = 1
                        x_red_use = x_red
                        stninfo_use = stninfo[gg, :]

                    tx_red = np.transpose(x_red_use)
                    twx_red = np.matmul(tx_red, w_pcp_red)

                    # calculate pcp
                    b = least_squares(x_red_use, y_prcp_red, twx_red)
                    pcpgg = np.dot(stninfo_use, b)
                    pcp_err_stn[gg, t] = pcpgg - y_prcp[gg]

            # tmean/trange processing
            if y_tmean[gg] > -100:
                # reduced matrices for tmean_stn/trange_stn
                nstn_temp = int(np.sum(near_stn_tempLoc[gg, :] > -1))
                if nstn_temp >= nearstn_min:
                    w_temp_red = np.zeros([nstn_temp, nstn_temp])
                    for i in range(nstn_temp):
                        w_temp_red[i, i] = near_stn_tempWeight[gg, i]  # eye matrix: stn weight in one-one lien
                    w_temp_1d = near_stn_tempWeight[gg, 0:nstn_temp]  # stn weight
                    w_temp_1d_loc = near_stn_tempLoc[gg, 0:nstn_temp]  # stn ID number/location
                    y_tmean_red = y_tmean[w_temp_1d_loc]  # transformed temp
                    y_trange_red = y_trange[w_temp_1d_loc]  # transformed temp
                    x_red_t = stninfo[w_temp_1d_loc, :]  # station lat/lon/ele/slope_ns/slope_we

                    ndata_t = np.sum(y_tmean_red > -100)
                    nodata_t = nstn_temp - ndata_t  # invalid temperature
                else:
                    x_red_t = 0  # not really necessary. just to stop warming from Pycharm
                    w_temp_red = 0  # not really necessary. just to stop warming from Pycharm
                    y_tmean_red = 0  # not really necessary. just to stop warming from Pycharm
                    y_trange_red = 0  # not really necessary. just to stop warming from Pycharm
                    ndata_t = 0
                    nodata_t = 0

                if ndata_t > 0:
                    stninfo_use = stninfo[gg, 0:4]
                    x_red_use = x_red_t[:, 0:4]  # do not use slope for temperature
                    tx_red = np.transpose(x_red_use)
                    twx_red = np.matmul(tx_red, w_temp_red)

                    b = least_squares(x_red_use, y_tmean_red, twx_red)
                    tmeangg = np.dot(stninfo_use, b)
                    tmean_err_stn[gg, t] = tmeangg - y_tmean[gg]

                    b = least_squares(x_red_use, y_trange_red, twx_red)
                    trangegg = np.dot(stninfo_use, b)
                    trange_err_stn[gg, t] = trangegg - y_trange[gg]

    return pcp_err_stn, tmean_err_stn, trange_err_stn


def regression(prcp_stn, tmean_stn, trange_stn, pcp_err_stn, tmean_err_stn, trange_err_stn, stninfo, gridinfo,
               mask, near_grid_prcpLoc, near_grid_prcpWeight, near_grid_tempLoc, near_grid_tempWeight,
               nearstn_min, nearstn_max, trans_exp, trans_mode):

    nstn, ntimes = np.shape(prcp_stn)
    nrows, ncols, nvars = np.shape(gridinfo)

    # 6.2 initialization
    tmp_weight_arr = np.eye(nearstn_max)
    y_max = -3 * np.ones([nrows, ncols, ntimes])
    pcp = -3 * np.ones([nrows, ncols, ntimes])
    pop = np.zeros([nrows, ncols, ntimes])
    pcp_err = np.zeros([nrows, ncols, ntimes])
    tmean = np.zeros([nrows, ncols, ntimes])
    trange = np.zeros([nrows, ncols, ntimes])
    tmean_err = np.zeros([nrows, ncols, ntimes])
    trange_err = np.zeros([nrows, ncols, ntimes])

    # start regression ...
    # loop through time steps
    for t in range(ntimes):
        print('Regression time step: ', t + 1, '---Total time steps: ', ntimes)
        # assign vectors of station alues for prcp_stn, temp, for current time step
        # transform prcp_stn to approximate normal distribution
        y_prcp = au.transform(prcp_stn[:, t], trans_exp, trans_mode)
        y_tmean = tmean_stn[:, t]
        y_trange = trange_stn[:, t]

        # loop through grids (row, col)
        for rr in range(nrows):
            for cc in range(ncols):
                if mask[rr, cc] != 1:
                    # the grid is outside mask extent
                    continue
                ########################################################################################################

                # 6.3 Precipitation estimation (pop and pcp)
                nstn_prcp = int(np.sum(near_grid_prcpLoc[rr, cc, :] > -1))
                if nstn_prcp < nearstn_min:
                    print('Precipitation regression: current time step, row, and col are', t, rr, cc)
                    sys.exit('Cannot find enough input stations for this grid cell')
                else:
                    # 6.3.1 reduced matrices for precipitation
                    w_pcp_red = np.zeros([nstn_prcp, nstn_prcp])
                    for i in range(nstn_prcp):
                        w_pcp_red[i, i] = near_grid_prcpWeight[rr, cc, i]  # eye matrix: stn weight in one-one lien
                    w_pcp_1d = near_grid_prcpWeight[rr, cc, 0:nstn_prcp]  # stn weight
                    w_pcp_1d_loc = near_grid_prcpLoc[rr, cc, 0:nstn_prcp]  # stn ID number/location
                    y_prcp_red = y_prcp[w_pcp_1d_loc]  # transformed prcp_stn
                    x_red = stninfo[w_pcp_1d_loc, :]  # station lat/lon/ele/slope_ns/slope_we

                    yp_red = np.zeros(nstn_prcp)  # pop: 0/1
                    yp_red[prcp_stn[w_pcp_1d_loc, t] > 0] = 1
                    ndata = np.sum(yp_red == 1)  # number of prcp_stn>0
                    nodata = np.sum(yp_red == 0)
                    y_max[rr, cc, t] = max(y_prcp_red)

                    # 6.3.2 estimate pop, pcp, pcp_err
                    # note: pcp_err is based on results from 6.1 and independent with gridded pop and pcp regression
                    if ndata == 0:
                        # nearby stations do not have positive prcp data (i.e., zero or missing)
                        pop[rr, cc, t] = 0
                        pcp[rr, cc, t] = y_prcp_red[0]  # corresponding to zero precipitation
                        pcp_err[rr, cc, t] = 0
                    else:
                        # decide if slope is to be used in regression
                        # tmp needs to be matmul(TX, X) where TX = TWX_red and X = X_red
                        mat_test = np.matmul(np.transpose(x_red), w_pcp_red)
                        tmp = np.matmul(mat_test, x_red)
                        vv = np.max(np.abs(tmp), axis=1)
                        if np.any(vv) == 0 or (abs(gridinfo[rr, cc, 4]) < 3.6 or abs(gridinfo[rr, cc, 5]) < 3.6):
                            slope_flag_pcp = 0
                            x_red_use = x_red[:, 0:4]  # do not use slope in regression
                            gridinfo_use = gridinfo[rr, cc, 0:4]
                        else:
                            slope_flag_pcp = 1
                            x_red_use = x_red
                            gridinfo_use = gridinfo[rr, cc, :]

                        tx_red = np.transpose(x_red_use)
                        twx_red = np.matmul(tx_red, w_pcp_red)

                        # calculate pop
                        if nodata == 0:
                            pop[rr, cc, t] = 1
                        else:
                            b = logistic_regression(x_red_use, twx_red, yp_red)
                            zb = - np.dot(gridinfo_use, b)
                            pop[rr, cc, t] = 1 / (1 + np.exp(zb))

                        # calculate pcp
                        b = least_squares(x_red_use, y_prcp_red, twx_red)
                        pcp[rr, cc, t] = np.dot(gridinfo_use, b)

                        # 6.4.3 estimate pcp error
                        err0 = pcp_err_stn[w_pcp_1d_loc, t]
                        pcp_err[rr, cc, t] = (np.sum((err0 ** 2) * w_pcp_1d) / np.sum(w_pcp_1d)) ** 0.5

                ############################################################################################################

                # 6.5 Temperature estimation (tmean and trange)
                # reduced matrices for tmean_stn/trange_stn
                nstn_temp = int(np.sum(near_grid_tempLoc[rr, cc, :] > -1))
                if nstn_temp < nearstn_min:
                    print('Temperature regression: current time step, row, and col are', t, rr, cc)
                    sys.exit('Cannot find enough input stations for this grid cell')
                else:
                    # 6.5.1 reduced matrices for precipitation
                    w_temp_red = np.zeros([nstn_temp, nstn_temp])
                    for i in range(nstn_temp):
                        w_temp_red[i, i] = near_grid_tempWeight[rr, cc, i]  # eye matrix: stn weight in one-one lien
                    w_temp_1d = near_grid_tempWeight[rr, cc, 0:nstn_temp]  # stn weight
                    w_temp_1d_loc = near_grid_tempLoc[rr, cc, 0:nstn_temp]  # stn ID number/location
                    y_tmean_red = y_tmean[w_temp_1d_loc]  # transformed temp
                    y_trange_red = y_trange[w_temp_1d_loc]  # transformed temp
                    x_red_t = stninfo[w_temp_1d_loc, :]  # station lat/lon/ele/slope_ns/slope_we

                    ndata_t = np.sum(y_tmean_red > -100)

                    if ndata_t == 0:
                        # This is not a problem for serially complete dataset (scd).
                        # But even if inputs have missing data, the way in Fortran-based GMET (simple filling) is not
                        # a good way. This problem should be solved when finding near stations.
                        print('Temperature regression: current time step, row, and col are', t, rr, cc)
                        sys.exit('Nearby stations do not have any valid temperature data')
                    else:
                        # 6.5.1 estimate tmean and its error
                        gridinfo_use = gridinfo[rr, cc, 0:4]
                        x_red_use = x_red_t[:, 0:4]  # do not use slope for temperature
                        tx_red = np.transpose(x_red_use)
                        twx_red = np.matmul(tx_red, w_temp_red)
                        b = least_squares(x_red_use, y_tmean_red, twx_red)
                        tmean[rr, cc, t] = np.dot(gridinfo_use, b)

                        # error estimation
                        err0 = tmean_err_stn[w_temp_1d_loc, t]
                        tmean_err[rr, cc, t] = (np.sum((err0 ** 2) * w_temp_1d) / np.sum(w_temp_1d)) ** 0.5

                        # 6.5.2 estimate trange and its error
                        b = least_squares(x_red_use, y_trange_red, twx_red)
                        trange[rr, cc, t] = np.dot(gridinfo_use, b)

                        # error estimation
                        err0 = trange_err_stn[w_temp_1d_loc, t]
                        trange_err[rr, cc, t] = (np.sum((err0 ** 2) * w_temp_1d) / np.sum(w_temp_1d)) ** 0.5

    return pop, pcp, tmean, trange, pcp_err, tmean_err, trange_err, y_max