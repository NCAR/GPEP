reapop_stn = np.nan * np.zeros([reanum, nstn, ntimes], dtype=np.float32)
timeall = np.hstack( (np.arange(1, 7000, 200), np.arange(7001, 7202, 50), np.arange(7401, 7600, 50), np.arange(7601, 14610, 200), 14610))
for tt in range(len(timeall) - 1):
    filet = path_pop + '/reapop_stn_' + str(timeall[tt]) + '-' + str(timeall[tt + 1]) + '.npz'
    d = np.load(filet)
    reapop_stn[:, :, (timeall[tt] - 1):timeall[tt + 1]] = d['reapop_stn'][:, :, (timeall[tt] - 1):timeall[tt + 1]]