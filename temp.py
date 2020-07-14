d = io.loadmat(outfile_stn)
data_stn = d['data_stn']
LLE = d['LLE']
date2 = d['date'][:]
del d