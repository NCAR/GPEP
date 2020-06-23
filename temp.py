# load near station information
datatemp = np.load(near_file_GMET)
if vars[v] == 'prcp' or vars[v] == 'pop':
    near_loc = datatemp['near_grid_prcpLoc']
    near_weight = datatemp['near_grid_prcpWeight']
    near_dist = datatemp['near_grid_prcpDist']
else:
    near_loc = datatemp['near_grid_tempLoc']
    near_weight = datatemp['near_grid_tempWeight']
    near_dist = datatemp['near_grid_tempDist']
near_loc = np.flipud(near_loc)
near_weight = np.flipud(near_weight)
near_dist = np.flipud(near_dist)
del datatemp