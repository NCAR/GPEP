import numpy as np
nstn = 27275
reapop_stn = np.nan * np.zeros([3, nstn, 14610], dtype=np.float32)
time = np.hstack([np.arange(1,14610,100),14610])
for i in range(147):
    print(i)
    filei = '/home/gut428/ReanalysisCorrMerge/pop/reapop_stn_' + str(time[i]) + '-' + str(time[i+1]) + '.npz'
    di = np.load(filei)
    reapop_stn[:, :, (time[i]-1):time[i+1]] = di['reapop_stn'][:, :, (time[i]-1):time[i+1]]
np.savez_compressed('/home/gut428/ReanalysisCorrMerge/pop/reanalysis_pop_stn.npz', reapop_stn=reapop_stn)