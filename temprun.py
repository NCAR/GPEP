import numpy as np
import sys
year = int(sys.argv[1])
path = '/home/gut428/OImerge_GWRLSBMA'


for m in range(12):
    file=path + '/oimerge_prcp' + str(year*100+m+1) + '.npz'
    d=np.load(file)
    oi_value = d['oi_value']
    oi_error = d['oi_error']
    latitude = d['latitude']
    longitude = d['longitude']
    oi_error = oi_error ** 0.5
    np.savez_compressed(file, oi_value=oi_value, oi_error=oi_error, latitude=latitude, longitude=longitude)

    file = path + '/oimerge_tmean' + str(year * 100 + m + 1) + '.npz'
    d = np.load(file)
    oi_value = d['oi_value']
    oi_error = d['oi_error']
    latitude = d['latitude']
    longitude = d['longitude']
    oi_error = oi_error ** 0.5
    np.savez_compressed(file, oi_value=oi_value, oi_error=oi_error, latitude=latitude, longitude=longitude)

    file = path + '/oimerge_trange' + str(year * 100 + m + 1) + '.npz'
    d = np.load(file)
    oi_value = d['oi_value']
    oi_error = d['oi_error']
    latitude = d['latitude']
    longitude = d['longitude']
    oi_error = oi_error ** 0.5
    np.savez_compressed(file, oi_value=oi_value, oi_error=oi_error, latitude=latitude, longitude=longitude)

    file = path + '/oimerge_prcp' + str(year * 100 + m + 1) + '_boxcox.npz'
    d = np.load(file)
    oi_value = d['oi_value']
    oi_error = d['oi_error']
    latitude = d['latitude']
    longitude = d['longitude']
    tranexp = d['tranexp']
    transmode = d['transmode']
    oi_error = oi_error ** 0.5
    np.savez_compressed(file, oi_value=oi_value, oi_error=oi_error, latitude=latitude, longitude=longitude, tranexp=tranexp, transmode=transmode)

    file = '/datastore/GLOBALWATER/CommonData/EMDNA_new/OImerge_pop/oimerge_pop' + str(year * 100 + m + 1) + '.npz'
    d = np.load(file)
    oi_value = d['oi_value']
    oi_error = d['oi_error']
    latitude = d['latitude']
    longitude = d['longitude']
    oi_error = oi_error ** 0.5
    file = '/home/gut428/OImerge_pop/oimerge_pop' + str(year * 100 + m + 1) + '.npz'
    np.savez_compressed(file, oi_value=oi_value, oi_error=oi_error, latitude=latitude, longitude=longitude)