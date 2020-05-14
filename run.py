import time
import numpy as np

data1=np.ones([30,14000])
weight=np.arange(30)

start = time.time()
data2=np.arange(14000)
for i in range(14000):
    datai = data1[:,i]
    data2[i] = np.sum(weight*datai)
end = time.time()
print(end-start)

start = time.time()
weight2 = np.tile(weight,(14000,1)).T
data3 = weight2*data1
data4 = np.nansum(data3,axis=0)
end = time.time()
print(end-start)

print(data4[1241],data2[1241])