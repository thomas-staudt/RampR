
import numpy as np
from scipy.integrate import quad

from . import core
from . import capi
from . import dist as distrib


#
# Helper function that allows to conduct experiments with a given xi(t)
# curve for the volume growth
#

def evolve_varying_xi( sim, xi, 
                       until_time = None,
                       until_vfrac = None,
                       initial_dt = None,
                       max_dt = 1.):
    
    if initial_dt: sim.target_dt = initial_dt

    nr_evaporated = 0

    if until_vfrac is not None:
        # Note: This works only approximately and might not stop at until_vfrac exactly
        while sim.vfrac < until_vfrac:
            sim.xi = xi(sim.time)
            if sim.xi > 0:
                until_time = sim.time + (until_vfrac-sim.vfrac) * sim.reference_volume / sim.xi
            else:
                until_time = sim.time + max_dt
            tmp = capi.lib.rampr_rkck_evolution_single(sim, min(until_time, sim.time + max_dt), 1)
            if tmp < 0: nr_evaporated = -1; break
            else: nr_evaporated += tmp

    else:
        while sim.time < until_time:
            sim.xi = xi(sim.time)
            tmp = capi.lib.rampr_rkck_evolution_single(sim, min(until_time, sim.time + max_dt), 1)
            if tmp < 0: nr_evaporated = -1; break
            else: nr_evaporated += tmp

    if nr_evaporated < 0:
        raise core.EvolutionError("%s" % sim.error_message)

    return nr_evaporated


#
# Convert radius rankings to probability density functions
#

def pdf(r):
    # This approximates the density -- should only be used for
    # visualization, not for further analysis
    grad = np.gradient(r, 1./len(r))
    density = -1. / grad
    # Norm it, as it is to be a prob.density
    density = density / (-1*np.trapz(density, r))

    x = np.flipud(r)
    y = np.flipud(density)

    # So that the density does not begin in the middle of the graph...
    y = np.hstack([[0.], y, [0.]])
    x = np.hstack([[x[0]], x, [x[-1]]])

    return x, y


#
# Convert radius rankings to cumulative density functions.
#

def cdf(r, append_value=0, prepend_value=None, include_evaporated=False):
    thr = 0. if not include_evaporated else 1. - float(nr_living)/float(r.size)

    # Get number of "living" particles
    nr_living = (r*r*r > capi.lib.EVAPORATION_TOLERANCE).sum()

    if append_value:  tmp = np.hstack( [ r[:nr_living], np.array([append_value]) ] )
    else:             tmp = r[:nr_living]

    if prepend_value: tmp = np.hstack( [ np.array([prepend_value]), tmp ] )

    x = np.flipud(tmp).copy()
    y = np.linspace(thr, 1, len(x))

    return x, y


#
# Get cumulative density of the LSW profile with order parameter p
# for given rho = r / rmax
#

def lswcdf(rho, p = np.inf):
    # Formula is taken from Niethammer, Pego (1999)
    if not isinstance(rho, (list, tuple, np.ndarray)):
        rho = np.array( [rho] )
    result = np.ones(rho.size)
    rho_ = rho[ rho < 1 ]
    if p == np.inf:
        result[rho < 1] = 1 - np.exp(-rho_ / (1-rho_)) /  \
               (np.power(1 - rho_, 5./3.) * np.power(1 + rho_/2., 4./3.))
    else:
        a  = 0.5 * (-1 + np.sqrt(3*(4*(1+1./p) - 1)))
        p1 = 3*a*a / ( (a-1) * (2*a+1) )
        p2 = 3*(a+1)*(a+1) / ((a+2) * (2*a+1))

        result[rho < 1] = 1 - np.power(1 - rho_, p) /  \
               ( np.power(1 - rho_/a, p1) * np.power(1 + rho_/(a+1), p2) )

    return result if result.size > 1 else result[0]


#
# Get cumulative density of the LSW profile with order parameter p
# for given rho = r / rmean
#

def lswcdf_mean(rho, p = np.inf):
    m = quad(lambda x: 1 - lswcdf(x, p), 0, 1)[0]
    return lswcdf(rho*m, p)


#
# Get density of LSW profile (p = infinity) with order parameter 
#

def lswpdf(rho):
    if isinstance(rho, (list, tuple, np.ndarray)):
        return np.array( [ lswpdf(rho_) for rho_ in rho ] )

    else:
        if rho >= 1.5 or rho < 0:
            return 0.

        a = 3. / (3. + rho)
        b = 3 / (3 - 2*rho)

        return 4./9. * (rho * rho) * a**(7./3.) * b**(11./3.) * np.exp(-b + 1)



_lsw_r_values   = np.linspace(0, 1.5, 10000)
_lsw_cdf_values = np.cumsum( lswpdf(_lsw_r_values) )
_lsw_cdf_values /= _lsw_cdf_values[-1]

def lswcdf_(rho):
    return np.interp(rho, _lsw_r_values, _lsw_cdf_values)


#
# Interpolate a radius distribution at the tip with a given order
#

def tip_interpolation(r, p):
    x = r[0] - r[1]
    if p != np.inf and p > 0:
        return r[0] - x / 2**(1./p)
    if p == np.inf:
        return r[0] - x / (1 + np.log(2) * x / r[0])
    else:
        raise ValueError('Interpolation order p must be positive or infinity')


#
# Calculate distances between radii distributions
#

def pdist(r, p):
    lsw_radii = distrib.lsw(1, len(r), p)
    return np.mean(np.abs(lsw_radii - r / np.mean(r)))


def tdist(r1, r2, n=1000):
    x1 = np.linspace(0, 1, len(r1))
    x2 = np.linspace(0, 1, len(r2))
    xtest = np.linspace(0, 1, n)
    r1test = np.interp(xtest, x1, r1)
    r2test = np.interp(xtest, x2, r2)

    return np.mean(np.abs(r1test - r2test))


def dist(r1, r2, n=1000):
    return tdist(r1 / np.mean(r1), r2 / np.mean(r2))

