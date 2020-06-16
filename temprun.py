import numpy as np
import calendar
from auxiliary_merge import m_DateList
nstn = 27275
pcp_reg_stn = np.nan * np.zeros([nstn, 14610], dtype=np.float32)
tmean_reg_stn = np.nan * np.zeros([nstn, 14610], dtype=np.float32)
trange_reg_stn = np.nan * np.zeros([nstn, 14610], dtype=np.float32)
pop_reg_stn = np.nan * np.zeros([nstn, 14610], dtype=np.float32)
date_list, date_number = m_DateList(1979, 2018, 'ByYear')

for y in range(1979,2019):
    for m in range(1,13):
        print(y,m)
        indym = (date_number['yyyy']==y) & (date_number['mm']==m)
        ed = calendar.monthrange(y,m)[1]
        str1 = str(y*10000+m*100+1)
        str2 = str(y*10000+m*100+ed)
        filei = '/home/gut428/GMET/PyGMETout_old/error_' + str1 + '-' + str2 + '.npz'
        di = np.load(filei)
        pcp_reg_stn[:,indym]=di['pcp_reg_stn']
        tmean_reg_stn[:, indym] = di['tmean_reg_stn']
        trange_reg_stn[:, indym] = di['trange_reg_stn']
        pop_reg_stn[:, indym] = di['pop_reg_stn']
np.savez_compressed('/home/gut428/GMET/PyGMETout_old/daily_regression_stn_old.npz', prcp=pcp_reg_stn,
                    tmean=tmean_reg_stn,trange=trange_reg_stn,pop=pop_reg_stn)