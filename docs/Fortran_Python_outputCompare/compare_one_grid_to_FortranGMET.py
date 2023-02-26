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

nearindexrc = tar_nearIndex[r, c, :]
stn_predictorrc = stn_predictor[nearindexrc, :]
nearWeightrc = tar_nearWeight[r, c, :]

# pcp
# fortran_order = np.array([113,692,616,613,664,663,629,597,638,659,641,658,585,637,683,108,119,1189,107,1142,99,106,105,668,103,660, 98,607,601,622,109,100,102, 97,101]) - 1
# tmean
# fortran_order = np.array([1046,1037,1059, 999,1142, 993,1010,1041, 989, 671,1013, 585,1003, 629, 622, 637,1049, 638, 619,1058, 664, 601,1024,1135, 658,1189,1184, 607, 691, 692, 660, 668,1052,1034,1185]) - 1
fortran_order = np.array([948,511,295,708,480,164,870,730,743,1197,983,855,303,321,817,1192,810,917,272,474,751,877,873,879,957,296,1150,228,794,535,883,227,280,887,789]) - 1

index_map = []
for i in fortran_order:
    index_map.append(np.where(nearindexrc==i)[0][0])

print(nearindexrc[index_map])
print(stn_predictorrc[index_map])
print(nearWeightrc[index_map])


print(ydata_tar)
print(xdata_near[index_map])
