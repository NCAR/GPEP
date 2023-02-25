# fortran
#,time step is,,1
#,grid number is, 66
import numpy  as np

# for debug mode
r = 134
c = 7

# target grid info is the same with Fortram
print(xaxis[c], yaxis[r]) # -120.65625, 34.53125
print('tar_predictor', tar_predictor[r, c, :])

# nearby station location index and predictors are the same (although the order of nearby stations are different)
# python distance is km  while fortran is n mi
# distance and weights are almost the same, with very small difference that can ignored, which comes from distance calculation difference

nearindexrc = nearIndex[r, c, :]
stn_predictorrc = stn_predictor[nearindexrc, :]
nearWeightrc = nearWeight[r, c, :]

nearDistancerc = ds_nearinfo.nearDistance_Grid_tmean.values[r, c, :]

# fortran_order = np.array([827,466,837,289,871,765,324,1155,878,1172,1205,838,399,437,770,727,1198,857,946,439,373,487,350,1204,1195,195,900,334,323,393,741,695,441,359,788]) - 1
fortran_order = np.array([113,692,616,613,664,663,629,597,638,659,641,658,585,637,683,108,119,1189,107,1142,99,106,105,668,103,660, 98,607,601,622,109,100,102, 97,101]) - 1

index_map = []
for i in fortran_order:
    index_map.append(np.where(nearindexrc==i)[0][0])

print(nearindexrc[index_map])
print(stn_predictorrc[index_map])
print(nearDistancerc[index_map])
print(nearDistancerc[index_map] * 0.539957) # km to n mi
print(nearWeightrc[index_map])


