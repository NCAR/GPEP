# perform locally weighted regression for all stations using leave-one-out


import numpy as np
import regression as reg
import os
import sys, time


x_red_use = np.array([ [1, 1.2, 1.4],
                       [1, 2.2, 1.42],
                       [1, 4.2, 3.33],
                       [1, 12, 3.5],
                       [1, 0.4, 2],
                       [1, 0.523, 1.2],
                       [1, 5.5, 2.3]])

w_pcp_red = np.zeros([7, 7])
w_vector = [0.1, 0.5, 0.23, 0.6, 0.23, 0.66, 0.9]
for i in range(7):
    w_pcp_red[i, i] = w_vector[i]  # eye matrix: stn weight in one-one line

y_prcp_red = np.array([0.2, 1.2, 0.4, 0.3, 0.55, 0.9, 0.44])

tx_red = np.transpose(x_red_use)
twx_red = np.matmul(tx_red, w_pcp_red)

t1 = time.time()
for i in range(100):
    b = reg.least_squares(x_red_use, y_prcp_red, twx_red)
t2 = time.time()
print(b)
print('time cost:', t2-t1)