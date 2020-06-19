# check the quality of station data based on regression estimates to ensure that some toxic stations are excluded
import numpy as np

filestn = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stndata_whole.npz'
filereg = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_regression/daily_regression_stn.npz'
file_nearstn = '/datastore/GLOBALWATER/CommonData/EMDNA_new/stn_regression/nearstn_catalog.npz'