
import numpy as np
from scipy.stats import gamma as sc_gamma
from scipy.stats import norm as sc_norm
from scipy.special import erfinv

from . import util

#
# Create a radius distribution with tip behavior of order p=1, 2
#

def linear(rmin, rmax, nr=1000):
    return np.linspace(rmin, rmax, nr)

def quadratic(rmin, rmax, nr=1000):
    y = np.linspace(0, 1, nr)
    delta = rmax - rmin
    return rmin + delta * (np.arccos(1 - 2*y) / np.pi)

#
# Normal distribution
#

def normal(mean, stddev, nr=1000, rmax=None):
    cdf = np.linspace(0, 1, nr+2)[1:-1]
    return sc_norm.ppf(cdf, loc=mean, scale=stddev)


#
# Gamma distribution
#

def gamma(shape, scale, nr=1000):
    cdf = np.linspace(0, 1, nr+2)[1:-1]
    return sc_gamma.ppf(cdf, shape, scale=scale)


#
# Take an arbitray density function f and get an (approximate) 
# radius vector out of it
#

def from_approx_pdf(f, r0, rmax, nr):
    import scipy.integrate as integrate
    a, _ = integrate.quad(f, r0, rmax)
    r = [r0]
    print(r)
    while r[-1] < rmax:
        r.append(a/(nr*f(r[-1])) + r[-1])

    return np.array(r)[::-1]



#
# Radius rankings corresponding to the family (f_p)
# of self-similar LSW-solutions
#

def lsw(mean, nr=1000, p=np.inf):
    y   = np.linspace(0, 1, nr)
    x   = np.linspace(0, 1.5, 1e6)
    cdf = util.lswcdf_mean(x, p=p)

    rranking = np.zeros(len(y))
    j = 0
    for i in range(0, len(y)):
        if j >= len(x): rranking[i] = mean * x[np.argmax(cdf)]
        while y[i] > cdf[j]: j += 1
        rranking[i] = mean * x[j] 

    return rranking[::-1]

