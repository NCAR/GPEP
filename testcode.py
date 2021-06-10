filein = OI_path + '/OI_' + str(year) + '.nc4'
fileout = pathy_out + '/OI_' + str(year) + '.nc4'
# print information
print('InFile:', filein)
print('OutFile:', fileout)
if os.path.isfile(fileout):
    print('Outfile exists. Continue ...')
    continue
# read raw EMDNA data
datain = read_EMDNA(filein)
# extract subset region
dataout = datain.sel(y=slice(latrange[1], latrange[0]), x=slice(lonrange[0], lonrange[1]))
# write data
dataout = dataout[['prcp', 'tmean', 'trange']]
dataout.to_netcdf(fileout, encoding=encoding)
del dataout, datain
