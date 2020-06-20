import numpy as np
from matplotlib import pyplot as plt

year = [1979, 2018]
month = [1, 12]
nrows = 800
ncols = 1300

outpath = '/home/gut428/figures/ERA5_day_ds'
inpath = '/datastore/GLOBALWATER/CommonData/EMDNA_new/ERA5_day_ds'
varname = 'data'

for m in range(month[0], month[1] + 1):
    print('month',m)
    nday = 0
    data = np.zeros([nrows, ncols])
    figm = outpath + '/month_' + str(m) + '.png'
    for y in range(year[0], year[1] + 1):
        file = inpath + '/ERA5_ds_tmean_' + str(y * 100 + m) + '.npz'
        d = np.load(file)
        dy = d[varname]
        data = data + np.sum(dy, axis=2)
        nday = nday + np.shape(dy)[2]

    data = data / nday
    plt.imshow(np.mean(data, axis=2), cmap='nipy_spectral')
    plt.colorbar(orientation='horizontal', shrink=0.6, pad=0.05)
    plt.clim([-40, 40])
    plt.axis('off')
    plt.savefig(figm, dpi=600)
