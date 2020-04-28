import numpy as np

def least_squares(x, y, tx):
    # In fortran version, ludcmp and lubksb are used to calcualte matrix inversion
    # call ludcmp(a, indx, d)
    # call lubksb(a, indx, b)

    # In Python version, numpy is used to calculate matrix inversion
    b = np.matmul(tx, y)
    a = np.matmul(tx, x)

    deta = np.linalg.det(a) # Compute the determinant of an array
    if deta ==0:
        print('Singular matrix')
        b[:] = 0
    else:
        ainv = np.linalg.inv(a)
        b = np.matmul(ainv, b)

    return b

def logistic_regression(x, tx, yp):
    # nstn: station number
    # nvars: station attributes (1, lat/lon/...), 1 is for regression
    nstn, nvars = np.shape(x)

    b = np.zeros(nvars) # regression coefficients (beta)
    p = np.zeros(nstn) # estimated pop (probability of precipitation)
    f = 0 # flag: continue or stop loops
    it = 0 # iteration times

    while f != 1:
        # check for divergence
        xb = -np.matmul(x, b)
        if np.any(xb > 50):
            f = 1
        else:
            p = 1 / (1 + np.exp(xb))

        # check for divergence
        if np.any(p > 0.9999):
            # logistic regression diverging
            f = 1
        else:
            v = np.zeros([nstn, nstn])  # diagonal variance matrix
            for i in range(nstn):
                v[i, i] = p[i] * (1 - p[i])
            xv = np.matmul(v, x)
            yn = yp - p # difference between station occurrence (0/1) and estimated pop: Bnew -Bold in Martyn and Slater 2006
            bn = least_squares (xv, yn, tx)

            # check: converging
            if np.any( np.abs(bn) > 1e-4 ):
                f = 0
            else:
                f = 1

            # check: iteration times
            if it > 20:
                f = 1
            else:
                f = 0

            b = b + bn # update coefficients

        it = it + 1

    return b

