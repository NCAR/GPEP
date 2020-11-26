# merge monthly outputs to yearly outputs

########################################################################################################################
# use xarray.merge. Don't know why this will need so large memory (>80G for a year)
# import xarray
# import glob
# import os, sys
#
# a = int(sys.argv[1])
# # b = int(sys.argv[2])
# year = [a, a]
#
# path = '/home/gut428/scratch/GMET/EMDNA_out/Estimate'
# # year = [1979, 1979]
#
# for y in range(year[0], year[1]+1):
#     for e in range(1, 2):
#         print('processing year/ens: ', y, e)
#         # merge 12 months
#         outfile = '{}/{}/EMDNA_{}.{:03d}.nc4'.format(path,y,y,e)
#         if not os.path.isfile(outfile):
#             filee = '{}/{}/*{:03d}.nc4'.format(path,y,e)
#             ds = xarray.merge([xarray.open_dataset(f) for f in glob.glob(filee)])
#             comp = dict(zlib=True, complevel=9)
#             encoding = {var: comp for var in ds.data_vars}
#             ds.to_netcdf(outfile, encoding=encoding)
#             for f in glob.glob(filee):
#                 os.remove(f)

########################################################################################################################
import netCDF4 as nc
import numpy as np
import os, datetime, sys

def dateymd(y):
    date_start = datetime.date(y, 1, 1)
    date_end = datetime.date(y, 12, 31)
    daynum = (date_end - date_start).days + 1
    date_ymd = np.zeros(daynum, dtype=int)
    dated = date_start
    for d in range(daynum):
        if d > 0:
            dated = dated + datetime.timedelta(days=1)
        date_ymd[d] = int(dated.strftime("%Y%m%d"))
    return date_ymd

def ncwrite(outfile, pcp, tmean, trange, lattar, lontar, yyyymmdd):
    nrows, ncols, ntimes = np.shape(pcp)
    ncfid = nc.Dataset(outfile, 'w', format='NETCDF4')
    ncfid.createDimension('y', nrows)
    ncfid.createDimension('x', ncols)
    ncfid.createDimension('time', ntimes)
    ncfid.createDimension('const', 1)
    varin = ncfid.createVariable('date', 'i8', ('time'), zlib=True, complevel=9)
    varin[:] = yyyymmdd
    varin.description = 'Date: yyyymmdd'
    varin = ncfid.createVariable('latitude', 'f4', ('y'), zlib=True, complevel=9)
    varin[:] = lattar
    varin.description = 'Center latitude of grids'
    varin = ncfid.createVariable('longitude', 'f4', ('x'), zlib=True, complevel=9)
    varin[:] = lontar
    varin.description = 'Center longitude of grids'
    varin = ncfid.createVariable('prcp', 'f4', ('time', 'x', 'y'), zlib=True, complevel=9, fill_value=-999)
    varin[:] = np.transpose(pcp, [2, 1, 0])
    varin.description = 'precipitation (mm/day)'
    varin = ncfid.createVariable('tmean', 'f4', ('time', 'x', 'y'), zlib=True, complevel=9, fill_value=-999)
    varin[:] = np.transpose(tmean, [2, 1, 0])
    varin.description = 'mean daily temperature (deg_C)'
    varin = ncfid.createVariable('trange', 'f4', ('time', 'x', 'y'), zlib=True, complevel=9, fill_value=-999)
    varin[:] = np.transpose(trange, [2, 1, 0])
    varin.description = 'daily temperature range (deg_C)'
    ncfid.close()


def ncwrite_oi(outfile, pcp, pop, tmean, trange, lattar, lontar, yyyymmdd):
    nrows, ncols, ntimes = np.shape(pcp)
    ncfid = nc.Dataset(outfile, 'w', format='NETCDF4')
    ncfid.createDimension('y', nrows)
    ncfid.createDimension('x', ncols)
    ncfid.createDimension('time', ntimes)
    ncfid.createDimension('const', 1)
    varin = ncfid.createVariable('date', 'i8', ('time'), zlib=True, complevel=9)
    varin[:] = yyyymmdd
    varin.description = 'Date: yyyymmdd'
    varin = ncfid.createVariable('latitude', 'f4', ('y'), zlib=True, complevel=9)
    varin[:] = lattar
    varin.description = 'Center latitude of grids'
    varin = ncfid.createVariable('longitude', 'f4', ('x'), zlib=True, complevel=9)
    varin[:] = lontar
    varin.description = 'Center longitude of grids'
    varin = ncfid.createVariable('prcp', 'f4', ('time', 'x', 'y'), zlib=True, complevel=9, fill_value=-999)
    varin[:] = np.transpose(pcp, [2, 1, 0])
    varin.description = 'precipitation (mm/day)'
    varin = ncfid.createVariable('pop', 'f4', ('time', 'x', 'y'), zlib=True, complevel=9)
    varin[:] = np.transpose(pop, [2, 1, 0])
    varin.description = 'Probability of precipitation (0-1)'
    varin = ncfid.createVariable('tmean', 'f4', ('time', 'x', 'y'), zlib=True, complevel=9, fill_value=-999)
    varin[:] = np.transpose(tmean, [2, 1, 0])
    varin.description = 'mean daily temperature (deg_C)'
    varin = ncfid.createVariable('trange', 'f4', ('time', 'x', 'y'), zlib=True, complevel=9, fill_value=-999)
    varin[:] = np.transpose(trange, [2, 1, 0])
    varin.description = 'daily temperature range (deg_C)'
    ncfid.close()

########################################################################################################################
# # probabilistic outputs: month to year
# year1 = int(sys.argv[1])
# year2 = int(sys.argv[2])
# ens1 = int(sys.argv[3])
# ens2 = int(sys.argv[4])
# year = [year1, year2]
# ens = [ens1, ens2]
# print('Setting year/ensemble', year, ens)
#
# # year = [1979, 1979]
# # ens = [1, 1]
# path = '/home/gut428/scratch/GMET/EMDNA_out/Estimate_spearman'
# lontar = np.arange(-180 + 0.05, -50, 0.1)
# lattar = np.arange(85 - 0.05, 5, -0.1)
#
# for y in range(year[0], year[1] + 1):
#     yyyymmdd = dateymd(y)
#
#     outfiley = '{}/{}/EMDNA_{}_mean.nc4'.format(path, y, y)
#     if os.path.isfile(outfiley):
#         continue
#     pcp_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)
#     tmean_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)
#     trange_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)
#
#     for e in range(ens[0], ens[1] + 1):
#         print('processing year/member', y, e)
#         # merge data from 12 months
#         outfilee = '{}/{}/EMDNA_{}.{:03d}.nc4'.format(path, y, y, e)
#         if os.path.isfile(outfilee):
#             fid = nc.Dataset(outfilee)
#             pcp_e = fid.variables['prcp'][:].data
#             pcp_e = np.transpose(pcp_e, [2, 1, 0])
#             tmean_e = fid.variables['tmean'][:].data
#             tmean_e = np.transpose(tmean_e, [2, 1, 0])
#             trange_e = fid.variables['trange'][:].data
#             trange_e = np.transpose(trange_e, [2, 1, 0])
#         else:
#             for m in range(1, 13):
#                 file = '{}/{}/EMDNA_{}.{:03d}.nc4'.format(path, y, y * 100 + m, e)
#                 ncfid = nc.Dataset(file)
#                 # read data
#                 pcpm = ncfid.variables['pcp'][:].data
#                 tmeanm = ncfid.variables['t_mean'][:].data
#                 trangem = ncfid.variables['t_range'][:].data
#                 ncfid.close()
#                 pcpm = np.transpose(pcpm, [1, 2, 0])
#                 pcpm = np.flipud(pcpm)
#                 tmeanm = np.transpose(tmeanm, [1, 2, 0])
#                 tmeanm = np.flipud(tmeanm)
#                 trangem = np.transpose(trangem, [1, 2, 0])
#                 trangem = np.flipud(trangem)
#                 # merge months
#                 if m == 1:
#                     pcp_e = pcpm.copy()
#                     tmean_e = tmeanm.copy()
#                     trange_e = trangem.copy()
#                 else:
#                     pcp_e = np.concatenate((pcp_e, pcpm), axis=2)
#                     tmean_e = np.concatenate((tmean_e, tmeanm), axis=2)
#                     trange_e = np.concatenate((trange_e, trangem), axis=2)
#                 del pcpm, tmeanm, trangem
#             # save outfile
#             ncwrite(outfilee, pcp_e, tmean_e, trange_e, lattar, lontar, yyyymmdd)
#             # for m in range(1, 13):
#             #     file = '{}/{}/EMDNA_{}.{:03d}.nc4'.format(path, y, y * 100 + m, e)
#             #     os.remove(file)
#
#         # yearly mean
#         pcp_e[pcp_e < -900] = np.nan
#         pcp_y = pcp_y + pcp_e
#         tmean_e[tmean_e < -900] = np.nan
#         tmean_y = tmean_y + tmean_e
#         trange_e[trange_e < -900] = np.nan
#         trange_y = trange_y + trange_e
#         del pcp_e, tmean_e, trange_e
#     # save outfile
#     pcp_y = pcp_y/100
#     tmean_y = tmean_y/100
#     trange_y = trange_y/100
#     pcp_y[np.isnan(pcp_y)] = -999
#     tmean_y[np.isnan(tmean_y)] = -999
#     trange_y[np.isnan(trange_y)] = -999
#     ncwrite(outfiley, pcp_y, tmean_y, trange_y, lattar, lontar, yyyymmdd)

########################################################################################################################
# # Optimal interpolation inputs: month to year
# year = [1979, 2018]
# path = '/home/gut428/scratch/GMET/GMET_OIinput'
# lontar = np.arange(-180 + 0.05, -50, 0.1)
# lattar = np.arange(85 - 0.05, 5, -0.1)
#
# for y in range(year[0], year[1] + 1):
#     print('processing year', y)
#     yyyymmdd = dateymd(y)
#     outfiley = '{}/OI_{}.nc4'.format(path, y)
#     if os.path.isfile(outfiley):
#         continue
#     pcp_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)
#     tmean_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)
#     trange_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)
#     pop_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)
#     for m in range(1, 13):
#         print('reading month', m)
#         file = '{}/reg_{}.nc'.format(path, y * 100 + m)
#         ncfid = nc.Dataset(file)
#         # read data
#         pcpm = ncfid.variables['pcp'][:].data
#         popm = ncfid.variables['pop'][:].data
#         tmeanm = ncfid.variables['tmean'][:].data
#         trangem = ncfid.variables['trange'][:].data
#         ncfid.close()
#         pcpm = np.transpose(pcpm, [1, 2, 0])
#         pcpm = np.flipud(pcpm)
#         popm = np.transpose(popm, [1, 2, 0])
#         popm = np.flipud(popm)
#         tmeanm = np.transpose(tmeanm, [1, 2, 0])
#         tmeanm = np.flipud(tmeanm)
#         trangem = np.transpose(trangem, [1, 2, 0])
#         trangem = np.flipud(trangem)
#         # transform
#         ind=popm >= 0
#         pcpm[ind] = (pcpm[ind]/3 + 1) ** 3
#         # merge months
#         if m == 1:
#             pcp_y = pcpm.copy()
#             tmean_y = tmeanm.copy()
#             trange_y = trangem.copy()
#         else:
#             pcp_y = np.concatenate((pcp_y, pcpm), axis=2)
#             tmean_y = np.concatenate((tmean_y, tmeanm), axis=2)
#             trange_y = np.concatenate((trange_y, trangem), axis=2)
#         del pcpm, tmeanm, trangem
#     # save outfile
#     print('save output file')
#     ncwrite_oi(outfiley, pcp_y, pop_y, tmean_y, trange_y, lattar, lontar, yyyymmdd)


########################################################################################################################
# probabilistic outputs: calculate spread
year1 = int(sys.argv[1])
year = [year1, year1]
print('Setting year', year)

path = '/home/gut428/projects/rpp-kshook/gut428/EMDNA_v1/Estimate'
lontar = np.arange(-180 + 0.05, -50, 0.1)
lattar = np.arange(85 - 0.05, 5, -0.1)

for y in range(year[0], year[1] + 1):
    yyyymmdd = dateymd(y)
    outfiley = 'EMDNA_{}_spread.nc4'.format(y)
    pcp_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)
    tmean_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)
    trange_y = np.zeros([800, 1300, len(yyyymmdd)], dtype=np.float32)

    infilemean = '{}/{}/EMDNA_{}_mean.nc4'.format(path, y, y)
    fid = nc.Dataset(infilemean)
    pcp_m = fid.variables['prcp'][:].data
    pcp_m = np.transpose(pcp_m, [2, 1, 0])
    tmean_m = fid.variables['tmean'][:].data
    tmean_m = np.transpose(tmean_m, [2, 1, 0])
    trange_m = fid.variables['trange'][:].data
    trange_m = np.transpose(trange_m, [2, 1, 0])
    fid.close()
    pcp_m[pcp_m < -900] = np.nan
    tmean_m[tmean_m < -900] = np.nan
    trange_m[trange_m < -900] = np.nan
    for e in range(1, 101):
        print('processing year/member', y, e)
        # merge data from 12 months
        infilee = '{}/{}/EMDNA_{}.{:03d}.nc4'.format(path, y, y, e)
        fid = nc.Dataset(infilee)
        pcp_e = fid.variables['prcp'][:].data
        pcp_e = np.transpose(pcp_e, [2, 1, 0])
        tmean_e = fid.variables['tmean'][:].data
        tmean_e = np.transpose(tmean_e, [2, 1, 0])
        trange_e = fid.variables['trange'][:].data
        trange_e = np.transpose(trange_e, [2, 1, 0])

        # yearly mean
        pcp_e[pcp_e < -900] = np.nan
        pcp_y = pcp_y + (pcp_e-pcp_m)**2
        tmean_e[tmean_e < -900] = np.nan
        tmean_y = tmean_y + (tmean_e-tmean_m)**2
        trange_e[trange_e < -900] = np.nan
        trange_y = trange_y + (trange_e-trange_m)**2
        del pcp_e, tmean_e, trange_e
    # save outfile
    pcp_y = (pcp_y/99)**0.5
    tmean_y = (tmean_y/99)**0.5
    trange_y = (trange_y/99)**0.5
    pcp_y[np.isnan(pcp_y)] = -999
    tmean_y[np.isnan(tmean_y)] = -999
    trange_y[np.isnan(trange_y)] = -999
    ncwrite(outfiley, pcp_y, tmean_y, trange_y, lattar, lontar, yyyymmdd)

