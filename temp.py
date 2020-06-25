for j in range(nstn):
    if np.isnan(stndata[j, 0]):
        continue
    csi[j, i] = cal_csi(stndata[j, :], reapop_stn[i, j, :], 0.5)
    bs[j, i] = cal_bs(stndata[j, :], reapop_stn[i, j, :])