# this single module can realize locally weighted linear regression
import numpy as np

def least_squares(x, y, tx):
    # In fortran version, ludcmp and lubksb are used to calcualte matrix inversion
    # call ludcmp(a, indx, d)
    # call lubksb(a, indx, b)
    # In Python version, numpy is used to calculate matrix inversion
    c = np.matmul(tx, y)
    a = np.matmul(tx, x)
    n = np.shape(a)[0]
    b = np.zeros(n)
    deta = np.linalg.det(a)  # Compute the determinant of an array
    if deta == 0:
        # print('Singular matrix')
        b[:] = 0
    else:
        ainv = np.linalg.inv(a)
        b = np.matmul(ainv, c)
    return b

def m_lwlr(nearinfo, weightnear, datanear, tarinfo):
    # # nearinfo: predictors from neighboring stations
    # # [station number, predictor number + 1] array with the first column being ones
    # nearinfo = np.zeros([nnum, npred+1])
    #
    # # weightnear: weight of neighboring stations
    # # [station number, station number] array with weights located in the diagonal
    # weightnear = np.zeros([nnum, nnum])
    # for i in range(nnum):
    #     weightnear[i, i] = 123
    #
    # # tarinfo:  predictors from target stations
    # # [predictor number + 1] vector with the first value being one
    # tarinfo = np.zeros(npred+1)
    #
    # # datanear: data from neighboring stations. [station number] vector
    # datanear = np.zeros(nnum)

    # start regression
    tx_red = np.transpose(nearinfo)
    twx_red = np.matmul(tx_red, weightnear)
    b = least_squares(nearinfo, datanear, twx_red)
    datatar = np.dot(tarinfo, b)
    return datatar
