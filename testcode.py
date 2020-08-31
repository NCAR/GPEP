    # prepare the mean monthly distribution of precipitation and temperature based ens and oi
    import numpy as np
    import netCDF4 as nc
    from scipy import io

    for year in range(2016, 2017):
        for month in range(6,7):
            ens_path = '/home/gut428/scratch/GMET/EMDNA_out/Estimate_spearman'
            ens_path = '{}/{}'.format(ens_path, year)
            # load oi data
            vars = ['pcp','t_mean','t_range']
            for i in range(3):
                datai = np.zeros([800, 1300, 100], dtype=np.float32)
                for e in range(1,101):
                    ens_file = '{}/EMDNA_{}.{:03d}.nc4'.format(ens_path, year * 100 + month, e)
                    d = nc.Dataset(ens_file)
                    di = d[vars[i]][:].data
                    datai[:, :, e-1] = di[0]
                    d.close()
            io.savemat('test.mat',{'data':datai}, do_compression=True)

