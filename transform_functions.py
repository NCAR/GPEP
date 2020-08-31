from numpy import log as ln
from numpy import exp
import numpy as np
import sys
from scipy.interpolate import interp1d


def bc(p, c):
    return (p ** c - 1) / c

def sfun(x, n, c1=None):
    # simon functions
    # x: precipitation data
    # n: function orders
    # c1: parameter 1
    # if n <=3, c1 is unnecessary
    if n > 3 and c1 == None:
        sys.exit('n>3 but c1 is abscent')
    if any(x == 0):
        sys.exit('input data contain zero values')

    if n == 1:
        return ln(x)
    if n == 2:
        return x + ln(x)
    if n == 3:
        return ln(x) + ln(x + 1)
    if n == 4:
        return x - x ** -c1
    if n == 5:
        return ln(x) - x ** -c1
    if n == 6:
        return ln(x) - (x + 1) ** -c1
    if n == 7:
        return (x + 1) ** -c1 + ln(x)
    if n == 8:
        return (x + 1) ** c1 + ln(x)
    if n == 9:
        return c1 * ln(x) + x
    if n == 10:
        return c1 * ln(x + 1) + ln(x)
    if n == 11:
        return x - (ln(x + 1)) ** -c1
    if n == 12:
        return ln(x) - (ln(x + 1)) ** -c1
    if n == 13:
        return ln(x + 1) - (ln(x + 1)) ** -c1
    if n == 14:
        return (ln(x + 1)) ** c1 + ln(x)
    if n == 15:
        return ln(x / c1 + 1) + ln(x)
    if n == 16:
        return x - ln(x ** -c1 + 1)
    if n == 17:
        return ln(x) - ln(x ** -c1 + 1)
    if n == 18:
        return ln(x + 1) - ln(x ** -c1 + 1)
    if n == 19:
        return ln(x ** c1 + 1) + ln(x)
    # Note: 20 is not from Simon's functions
    if n == 20:
        return bc(x, c1)

def sfun_inv(x, n, c1=None):
    if n > 3 and c1 == None:
        sys.exit('n>3 but c1 is abscent')
    if any(x == 0):
        sys.exit('input data contain zero values')

    if n == 1:
        return exp(x)
    elif n <= 19:
        return sfun_invint(x, n, c1)
    if n == 20:
        return (x * c1 + 1) ** (1/c1)

def sfun_invint(x, n, c1):
    y0 = np.arange(0.01, 500, 0.1)
    x0 = sfun(y0, n, c1)
    f = interp1d(x0, y0)
    x[x < x0[0]] = x0[0]
    x[x > x0[-1]] = x0[-1]
    y = f(x)
    return y

