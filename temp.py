transmode = 'box-cox'
tranexp = 3
fileoi = path_oi + '/oimerge_prcp' + str(year * 100 + month) + '.npz'
datatemp = np.load(fileoi)
prcp = datatemp['oi_value']  # value in normal space
for d in range(ntimes):
    prcp[:, :, d] = prcp[:, :, d] * corr_ratio[:, :, month-1]
prcp = transform(prcp, tranexp, transmode)

fileoi = path_oi + '/oimerge_prcp' + str(year * 100 + month) + '_boxcox.npz'
datatemp = np.load(fileoi)
prcp_raw = datatemp['oi_value']  # value in normal space
prcp_err = datatemp['oi_error']
temp1 = (prcp + 3)
temp2 = (prcp_raw + 3)
ratio_err = temp1 / temp2  # mean value change
ratio_err[temp2 == 0] = 1
ratio_err[ratio_err > 1.5] = 1.5
ratio_err[ratio_err < 0.5] = 0.5
prcp_err = prcp_err * ratio_err