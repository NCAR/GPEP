# based on the OI merged precipitation and error estimates (natural space), apply transformation to approximate normal distribution
# this is temporarily named as 8.6, but is coded after s10_month2year. After everything is finished, those need to be renumbered

# results show this step cannot achieve better results than Box-Cox with a fixed parameter
# s8.6 and s8.7 are not used in the final product

import numpy as np
from scipy import io
import os, sys, time
import netCDF4 as nc
from transform_functions import sfun
from scipy import optimize
import lmoments3 as lm
import multiprocessing as mp
from tqdm import tqdm

def of_fun(param, p, num):
    pt = sfun(p, num, param)
    # lmoms = lmom.samlmu(pt, nmom=4)
    lmoms = lm.lmom_ratios(pt, nmom=4)
    l_skewness = lmoms[2]
    l_kurtosis = lmoms[3]
    of = (l_skewness) ** 2 + (l_kurtosis - 0.1226) ** 2
    return of

def of_solve_para(i, pi, tn1, tn2):
    out = np.zeros([tn2-tn1+1,2])
    # for t in range(3, 20):
    for t in range(tn1-1, tn2):
        paramrc = optimize.fmin(of_fun, 0.1, (pi, t + 1), disp=False)
        out[t-tn1+1, 0] = paramrc
        out[t-tn1+1, 1] = of_fun(paramrc, pi, t + 1)
    return (i, out)

########################################################################################################################

# time periods and methods
vars = sys.argv[1]  # prcp only
vars = [vars]
month = int(sys.argv[2])
transn1 = int(sys.argv[3])
transn2 = int(sys.argv[4])

# vars = 'prcp'
# vars = [vars]
# month = 1
# transn1 = 10
# transn2 = 10

print(vars, month, transn1, transn2)

########################################################################################################################

# ### Local Mac settings
# # input files/paths
# FileGridInfo = '/Users/localuser/Research/EMDNA/basicinfo/gridinfo_whole.nc'
# path_bac = '/Users/localuser/Research/EMDNA/merge' # data that will be used as background
# path_obs = '/Users/localuser/Research/EMDNA/regression' # data that will be used as observation
# near_file_GMET = '/Users/localuser/Research/EMDNA/regression/weight_nearstn.npz' # near station of stations/grids
# file_mask = './DEM/NA_DEM_010deg_trim.mat'
# FileStnInfo = '/Users/localuser/Research/EMDNA/basicinfo/stnlist_whole.txt'
# gmet_stndatafile = '/Users/localuser/Research/EMDNA/stndata_whole.npz'
#
# # output files/paths (can also be used as inputs once generated)
# path_oimerge = '/Users/localuser/Research/EMDNA/oimerge'
# ### Local Mac settings


# ### Plato settings
# # input files/paths
# FileGridInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/gridinfo_whole.nc'
# FileStnInfo = '/datastore/GLOBALWATER/CommonData/EMDNA_new/StnGridInfo/stnlist_whole.txt'
# gmet_stndatafile = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_aftercheck.npz'
# path_bac = '/datastore/GLOBALWATER/CommonData/EMDNA_new/ReanalysisCorrMerge/pop'
# path_obs = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck'
# near_file_GMET = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_reg_aftercheck/nearstn_catalog.npz'
# file_mask = '/datastore/GLOBALWATER/CommonData/EMDNA_new/DEM/NA_DEM_010deg_trim.mat'


### Graham settings
# input files/paths
FileGridInfo = '/home/gut428/projects/rpp-kshook/gut428/EMDNA_new/StnGridInfo/gridinfo_whole.nc'
gmet_stndatafile = '/home/gut428/projects/rpp-kshook/gut428/EMDNA_new/stndata_aftercheck.npz'
file_mask = '/home/gut428/projects/rpp-kshook/gut428/EMDNA_new/DEM/NA_DEM_010deg_trim.mat'
path_oimerge = '/home/gut428/projects/rpp-kshook/gut428/EMDNA_new/OImerge_GWRLSBMA'
outpath = '/home/gut428/scratch/OIprcp_transform'

# settings
step = 1 # all data in one 10X10 window are assembled
tnum = transn2 - transn1 + 1  # number of transformation. Coded as a constant parameter in of_solve_para
ncpus = mp.cpu_count() # for parallel
print('CPU number is ', ncpus)


########################################################################################################################

# basic processing
mask = io.loadmat(file_mask)
mask = mask['DEM']
mask[~np.isnan(mask)] = 1  # 1: valid pixels
nrows, ncols = np.shape(mask)

# meshed lat/lon of the target region
ncfid = nc.Dataset(FileGridInfo)
lattarm = ncfid.variables['latitude'][:].data
lattarm = np.flipud(lattarm)
lontarm = ncfid.variables['longitude'][:].data
ncfid.close()
lontar = lontarm[0, :]
lattar = lattarm[:, 0]

# load observations for all stations
datatemp = np.load(gmet_stndatafile)
date_ymd = datatemp['date_ymd']
del datatemp
date_yyyy = (date_ymd / 10000).astype(int)
date_mm = (np.mod(date_ymd, 10000) / 100).astype(int)

########################################################################################################################
# OI-merging at grid scale
for v in range(len(vars)):
    for m in range(month - 1, month):
        print('month', m + 1)
        outfile = outpath + '/prcp_trans_param_' + str(m+1) + '_TN' + str(transn1) + '-' + str(transn2) + '.mat'
        if os.path.isfile(outfile):
            print('outfile exists')
            continue
        indm = (date_mm == m + 1)
        nday = sum(indm)
        datem = date_yyyy[indm]

        # load all data for this month
        oi_value = np.nan * np.zeros([nrows, ncols, nday], dtype=np.float32)
        # oi_error = np.nan * np.zeros([nrows, ncols, nday], dtype=np.float32)
        for y in range(1979, 2019):
            print('year', y)
            fileoi_ym = path_oimerge + '/oimerge_' + vars[v] + str(y * 100 + m + 1) + '.npz'
            d = np.load(fileoi_ym)
            indym1 = datem == y
            oi_value[:, :, indym1] = d['oi_value']
            # oi_error[:, :, indym1] = d['oi_error']

        # calculate transformation parameters
        # normal version: to slow. about 2s for one r/c and all t
        param = np.nan * np.zeros([nrows, ncols, tnum, 2], dtype=np.float32)
        totnum = np.sum(~np.isnan(mask))
        flag = 1
        for r in range(nrows):
            for c in range(ncols):
                if not np.isnan(mask[r, c]):
                    if np.mod(flag, 100) == 0:
                        print('processing num/totnum/ratio:',flag, totnum, flag/totnum)
                    pi = oi_value[r, c, :]
                    pi = pi[pi > 0]
                    if len(pi) > 100:
                        for t in range(transn1 - 1, transn2):
                            paramrc = optimize.fmin(of_fun, 0.1, (pi, t + 1), disp=False)
                            param[r, c, t - transn1 + 1, 0] = paramrc
                            param[r, c, t - transn1 + 1, 1] = of_fun(paramrc, pi, t + 1)
                    flag = flag + 1
        # parallel
        # nrows2 = int(nrows/step)
        # ncols2 = int(ncols/step)
        # index2 = []
        # data2 = []
        # for r in range(nrows2):
        #     for c in range(ncols2):
        #         drc = oi_value[r*step:(r+1)*step, c*step:(c+1)*step, :]
        #         drc = drc.flatten()
        #         drc = drc[drc>0]
        #         if len(drc) > 100:
        #             data2.append(drc)
        #             index2.append([r, c])
        # del oi_value
        #
        # ns=np.shape(data2)[0]
        # print('start estimating parameters')
        # pool = mp.Pool(ncpus)
        # results = []
        # # call apply_async() without callback
        # pbar = tqdm(total=ns)
        # result_objects = [pool.apply_async(of_solve_para, args=(i, row, transn1, transn2), callback=lambda _: pbar.update(1))
        #                   for i, row in enumerate(data2)]
        # # result_objects is a list of pool.ApplyResult objects
        # results = [r.get()[1] for r in result_objects]
        #
        # param = np.nan * np.zeros([nrows, ncols, transn2-transn1+1, 2], dtype=np.float32)
        # for i in range(ns):
        #     index2i = index2[i]
        #     for r in range(index2i[0]*step, index2i[0]*step+step):
        #         for c in range(index2i[1]*step, index2i[1]*step+step):
        #             param[r, c, :, :] = results[i]
        # print('end estimating parameters')
        #
        # for t in range(tnum):
        #     ti = param[:, :, t, 0].copy()
        #     ti[np.isnan(ti)] = np.nanmean(ti)
        #     ti[np.isnan(mask)] = np.nan
        #     param[:, :, t, 0] = ti
        #
        # pool.close()
        # pool.join()
        io.savemat(outfile,{'param':param,'lontar':lontar,'lattar':lattar}, do_compression=True)

