var = varinfo[0]
file = '{}/{}{}.nc'.format(inpath, var['folder'], y * 100 + m + 1)
fid = nc.Dataset(file)
