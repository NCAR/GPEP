import numpy as np
from matplotlib import pyplot as plt
year = [1979, 2018]
month = [1, 12]
nrows = 800
ncols = 1300
vars=['prcp','tmean','trange']
clims=[[0,10], [-40,40], [0,20]]

for v in range(len(vars)):
    var = vars[v]
    varname = 'oi_value'
    outpath = '/home/gut428/figures/OImerge_GWRLSBMA'
    inpath = '/datastore/GLOBALWATER/CommonData/EMDNA_new/OImerge_GWRLSBMA'
    # monthly
    for y in range(2018, 2019):
        print('year', y)
        for m in range(month[0], month[1] + 1):
            file = inpath + '/oimerge_' + var + str(y * 100 + m) + '.npz'
            figm = outpath + '/' + var + '_' + str(y * 100 + m) + '.png'
            d = np.load(file)
            data = d[varname]
            data = np.mean(data, axis=2)
            plt.imshow(data, cmap='nipy_spectral')
            plt.colorbar(orientation='horizontal', shrink=0.6, pad=0.05)
            plt.clim(clims[v])
            plt.title('Date:'+str(y * 100 + m))
            plt.axis('off')
            plt.savefig(figm, dpi=600)
            plt.close()
    # monthly climatology
    for m in range(month[0], month[1] + 1):
        print('month',m)
        nday = 0
        data = np.zeros([nrows, ncols])
        figm = outpath + '/' + var + '_' + str(m) + '.png'
        for y in range(year[0], year[1] + 1):
            file = inpath + '/oimerge_' + var + str(y * 100 + m) + '.npz'
            d = np.load(file)
            dy = d[varname]
            data = data + np.sum(dy, axis=2)
            nday = nday + np.shape(dy)[2]
        data = data / nday
        plt.imshow(data, cmap='nipy_spectral')
        plt.colorbar(orientation='horizontal', shrink=0.6, pad=0.05)
        plt.clim(clims[v])
        plt.axis('off')
        plt.savefig(figm, dpi=600)
        plt.close()


for v in range(len(vars)):
    var = vars[v]
    varname = 'oi_value'
    outpath = '/home/gut428/figures/OImerge_GWRLSBMA'
    inpath = '/datastore/GLOBALWATER/CommonData/EMDNA_new/OImerge_GWRLSBMA'
    # monthly
    for y in range(2018, 2019):
        print('year', y)
        for m in range(month[0], month[1] + 1):
            file = inpath + '/oimerge_' + var + str(y * 100 + m) + '.npz'
            figm = outpath + '/' + var + '_' + str(y * 10000 + m * 100 +1) + '.png'
            d = np.load(file)
            data = d[varname]
            data = data[:, :, 0]
            plt.imshow(data, cmap='nipy_spectral')
            plt.colorbar(orientation='horizontal', shrink=0.6, pad=0.05)
            plt.clim(clims[v])
            plt.title('Date:'+str(y * 10000 + m * 100 + 1))
            plt.axis('off')
            plt.savefig(figm, dpi=600)
            plt.close()