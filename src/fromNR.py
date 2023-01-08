# "Numerical Recipes" ludcmp() C code on page 46 translated into Python
# I make it as simple as possible, disregard efficiency.
# Here a is a list of list, n is integer (size of the matrix)
# index is a list, and d is also a list of size 1
# Python list index starts from 0.  So matrix index is from 0 to n-1.
#
import math

def ludcmp(a, n, indx, d):
    d[0] = 1.0
    # looking for the largest a in each row and store it in vv as inverse
    # We need a new list same size as indx, for this we use .copy()
    vv = indx.copy()
    for i in range(0, n):
        big = 0.0
        for j in range(0, n):
            temp = math.fabs(a[i][j])
            if (temp > big):
                big = temp
        vv[i] = 1.0 / big
    #
    # run Crout's algorithm
    for j in range(0, n):
        # top half & bottom part are combined
        # but the upper limit l for k sum is different
        big = 0.0
        for i in range(0, n):
            if (i < j):
                l = i
            else:
                l = j
            sum = a[i][j]
            for k in range(0, l):
                sum -= a[i][k] * a[k][j]
            # end for k
            a[i][j] = sum
            # for bottom half, we keep track which row is larger
            if (i >= j):
                dum = vv[i] * math.fabs(sum)
                if (dum >= big):
                    big = dum
                    imax = i
            # end if (i>= ...)
        # end for i
        # pivoting part, swap row j with row imax, a[j] is a whole row
        if (j != imax):
            dum = a[imax]
            a[imax] = a[j]
            a[j] = dum
            d[0] = - d[0]
            vv[imax] = vv[j]
        # end if (j != ...)
        # divide by the beta diagonal value
        indx[j] = imax
        dum = 1.0 / a[j][j]
        for i in range(j + 1, n):
            a[i][j] *= dum
        # end for i
    # end for j


# end of def ludcmp

# We do backward substitution in lubksb() take the row swapped LU decomposed
# a, size n, and swapping indx, and b vector as input.  The output is
# in b after calling.
def lubksb(a, n, indx, b):
    ii = -1
    # forward
    for i in range(0, n):
        ip = indx[i]
        sum = b[ip]
        b[ip] = b[i]
        if (ii != -1):
            for j in range(ii, i):
                sum -= a[i][j] * b[j]
        elif (sum != 0):
            ii = i
        b[i] = sum
    # bote alpha_{ii} is 1 above
    #  backward
    for i in range(n - 1, -1, -1):
        sum = b[i]
        for j in range(i + 1, n):
            sum -= a[i][j] * b[j]
        b[i] = sum / a[i][i]


# end lubksb()


# unfortunately a is destroyed (become swapped LU)
def linearsolver(a, n, b):
    indx = list(range(n))
    d = [1]
    ludcmp(a, n, indx, d)
    x = b.copy()
    lubksb(a, n, indx, x)
    # print("x=", x)
    return x

if __name__ == "__main__":
    # end linearsolver
    a = [[1, 3, 3, -5], [2, -4, 7, -1], [7, 0.5, 3, -6], [9, -2, 3, 8]]
    n = 4
    b = [0, 2, 3, -10]
    # using ludcmp and lubksb
    linearsolver(a, n, b)

    # just using numpy
    from numpy.linalg import inv
    import numpy as np
    a = np.array([[1, 3, 3, -5], [2, -4, 7, -1], [7, 0.5, 3, -6], [9, -2, 3, 8]])
    b = np.array([0, 2, 3, -10])
    ainv = inv(a)
    x2 = np.matmul(ainv, b)
    x2