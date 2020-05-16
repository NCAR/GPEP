import numpy as np

def OImerge(tar_err_b, near_err_b, near_err_o):
    # tar/near: target and nearby stations/grids
    # o/b: observation/background
    # err: error (o-t or b-t where t is truth)
    # calculate weight (W) using: W(Cb + Co)=Cb0
    ntimes, nnum = np.shape(near_err_b)
    # covariance matrix of errors
    Cb = np.cov(near_err_b)
    Co = np.cov(near_err_o)
    Cb0 = np.zeros(nnum)
    for i in range(nnum):
        Cb0[i] = cov1d(tar_err_b, near_err_b[:,i])
    cbot = np.linalg.inv(Cb+Co)
    W = np.dot(Cb0, cbot)
    return W

def cov1d(a, b):
    if len(a) != len(b):
        return
    a_mean = np.mean(a)
    b_mean = np.mean(b)

    sum = 0

    for i in range(0, len(a)):
        sum += ((a[i] - a_mean) * (b[i] - b_mean))

    return sum/(len(a)-1)