
import numpy as np

from . import capi


# Error that may occur during calling evolution functions from the c-api
class EvolutionError(Exception): 
    pass

# Error that may occur during calling functions from the c-api that change
# the droplet distribution
class DistributionError(Exception):
    pass


# Need to inherit from object to make @properties work correctly
class Simulation(object):

    #
    # Constructor and destructor methods
    # A destructor is needed to prevent memory leaks
    #

    def __init__(self, radii, k, t0 = 0.):

        # Make sure that we can pass the radii array to the c init function
        # -- it has to be contigous in memory
        self.initial_radii = np.ascontiguousarray(radii)


        # Now, call the init function and obtain
        # the void ptr 'csimptr'. 
        self.csimptr = capi.lib.rampr_init_simulation( self.initial_radii, 
                                                        self.initial_radii.size, 
                                                        k )

        # Check if memory allocation went fine -- if it didnt, self.csimptr
        # is NULL, thus None in ctypes
        if self.csimptr == None:
            raise MemoryError( "Library '%s': Could not allocate memory "   + 
                               "in rampr_init_simulation." % capi.lib_name )


        # The _as_parameter_ member is for utility: We can then pass the
        # simulation instance sim itself to cytpes wrapped function,
        # instead of passing sim.csimptr
        self._as_parameter_ = self.csimptr

        
        # Set the initial time value via @setters
        self.time = t0


    def __del__(self):
        # free the memory allocated by librampr.so
        capi.lib.rampr_free_simulation(self)



    #
    # The core evolve functions: let the droplets evolve (for different
    # routines).
    # 

    def evolve_euler( self, until_vfrac = None, 
                      until_time = None, 
                      until_nr   = 0,
                      steps = 1, initial_dt = None ):
        # Only one of the three of until_vfrac, until_time, steps should be
        # specified when calling the function. Elsewise, until_vfrac will
        # preceed until_time will preceed steps.

        # If the evolution shall start with a given initial timestep dt
        # (and not the one stored in the target_dt member of the
        # simulation), update the target_dt
        if initial_dt: self.target_dt = initial_dt

        # Count the number of evaporated particles
        nr_evaporated = 0

        # Convert vfrac to time point
        if until_vfrac: 
            if self.k == 1: raise EvolutionError("Can't use until_vfrac > 1 for k = 1.")
            until_time = (until_vfrac - 1) * self.reference_volume / self.xi
        
        if until_nr and not until_time:
            raise EvolutionError("Must also specify until_time when using until_nr")

        # Evolve until the right time, if it is given
        if until_time:
            nr_evaporated = int( capi.lib.rampr_euler_evolution_until(self, until_time, until_nr, 0) )

        # Default: make one step. Can specify more steps, however
        else:
            for step in range(steps):
                tmp = int( capi.lib.rampr_euler_evolution_single(self, -1., 0) )
                if tmp < 0:
                    nr_evaporated = -1
                    break
                else:
                    nr_evaporated += tmp


        if nr_evaporated < 0: 
            raise EvolutionError("%s" % self.error_message)

        return nr_evaporated



    def evolve_rkck( self, until_vfrac = None, 
                     until_time = None, 
                     until_nr   = 0,
                     steps = 1, initial_dt = None, 
                     error_tolerance = None,
                     adapted = True
                   ):

        if initial_dt:      self.target_dt       = initial_dt
        if error_tolerance: self.error_tolerance = error_tolerance

        # Count the number of evaporated particles
        nr_evaporated = 0

        # Convert vfrac to time point
        if until_vfrac: 
            if self.k == 1: raise EvolutionError("Can't use until_vfrac > 1 for k = 1.")
            until_time = (until_vfrac - 1) * self.reference_volume / self.xi
        

        if until_nr and not until_time:
            raise EvolutionError("Must also specify until_time/until_vfrac when using until_nr")

        # Evolve until the right time, if it is given
        if until_time:
            nr_evaporated = int( capi.lib.rampr_rkck_evolution_until(self, until_time, until_nr, adapted) )
        # Default: make one step. Can specify more steps, however
        else:
            for step in range(steps):
                tmp = int( capi.lib.rampr_rkck_evolution_single(self, -1., adapted) )
                if tmp < 0:
                    nr_evaporated = -1
                    break
                else:
                    nr_evaporated += tmp


        if nr_evaporated < 0: 
            raise EvolutionError("%s" % self.error_message)

        return nr_evaporated

    def revive(self, radii, preserve_vfrac=True, preserve_k=True):
        radii = np.ascontiguousarray(radii)
        if not capi.lib.rampr_revive_droplets( self, radii, radii.size, preserve_vfrac, preserve_k ):
            raise DistributionError("%s" % self.error_message)

    def add(self, radii, preserve_vfrac=True, preserve_k=True, ignore_dead=False):
        radii = np.ascontiguousarray(radii)
        if not capi.lib.rampr_add_droplets( self, radii, radii.size, preserve_vfrac, preserve_k, ignore_dead ):
            raise DistributionError("%s" % self.error_message)



    # Overview string over the current status of the distribution
    def status(self, comment_char=None):
        stat = ( "RampR simulation with %d droplets in total at time t = %.3f\n\n"
                 "  k     = %10.3f    rmax  = %6.3f    rhomax  = %.3f\n"
                 "  thr   = %10.3f    rmin  = %6.3f    rhomin  = %.3f\n"
                 "  V/V0  = %10.3f    rmean = %6.3f    rhomean = %.3f\n"
                 "                        rcrit = %6.3f    rhocrit = %.3f\n"
               ) % (self.nr_total, self.time, 
                    self.k,          self.rmax,  self.rhomax, 
                    self.thr,        self.rmin,  self.rhomin,
                    self.vfrac,      self.rmean, self.rhomean, 
                                     self.rcrit, self.rhocrit)

        # If the status is to be commented, place comment_char at the
        # beginning of stat and after each newline
        if comment_char != None: 
            stat = comment_char + stat.replace("\n", "\n" + comment_char + " ")
        return stat

    # How shall the class be represented in the REPL, or when printing?
    def __repr__(self):
        return self.status()


    # Properties of the Simulation class. More convenient, automatically
    # call the underlying c-routines
    @property
    def nr_dead(self):
        return int(capi.lib.rampr_get_nr_dead(self))

    @property
    def nr_living(self):
        return int(capi.lib.rampr_get_nr_living(self))

    @property
    def nr_total(self):
        return int(capi.lib.rampr_get_nr_total(self))

    @property
    def time(self):
        return capi.lib.rampr_get_time(self)

    @time.setter
    def time(self, val):
        capi.lib.rampr_set_time(self, val)


    @property
    def radii(self):
        return np.ctypeslib.as_array(capi.lib.rampr_get_radii(self), (self.nr_total,))

    @property
    def r(self):
        return self.radii

    @property
    def sum_radii(self):
        return capi.lib.rampr_get_sum_radii(self)

    @property
    def rho(self):
        return self.radii/self.radii[0]

    @property
    def volumes(self):
        return np.ctypeslib.as_array(capi.lib.rampr_get_volumes(self), (self.nr_total,))

    @property
    def v(self):
        return self.volumes

    @property
    def sum_volumes(self):
        return capi.lib.rampr_get_sum_volumes(self)

    @property
    def reference_volume(self):
        return capi.lib.rampr_get_reference_volume(self)

    @reference_volume.setter
    def reference_volume(self, val):
        return capi.lib.rampr_set_reference_volume(self, val)


    @property
    def xi(self):
        return capi.lib.rampr_get_volume_growth_rate(self)

    @xi.setter
    def xi(self, val):
        capi.lib.rampr_set_volume_growth_rate(self, val)


    @property
    def k(self):
        return capi.lib.rampr_get_k_value(self)

    @k.setter
    def k(self, val):
        capi.lib.rampr_set_k_value(self, val)

    @property
    def last_dt(self):
        return capi.lib.rampr_get_last_dt(self)

    @property
    def target_dt(self):
        return capi.lib.rampr_get_target_dt(self)

    @target_dt.setter
    def target_dt(self, val):
        return capi.lib.rampr_set_target_dt(self, val)

    @property
    def error_tolerance(self):
        return capi.lib.rampr_get_error_tolerance(self)

    @error_tolerance.setter
    def error_tolerance(self, val):
        return capi.lib.rampr_set_error_tolerance(self, val)

    @property
    def nr_failed_adaptions(self):
        return int( capi.lib.rampr_get_nr_failed_adaptions(self) )

    @property
    def nr_steps(self):
        return int( capi.lib.rampr_get_nr_steps(self) )

    @property
    def rmin(self):
        return self.r[self.nr_living - 1]

    @property
    def rmax(self):
        return self.r[0]

    @property
    def rmean(self):
        return self.sum_radii / self.nr_living

    @property
    def vmin(self):
        return self.v[self.nr_living - 1]

    @property
    def vmax(self):
        return self.v[0]

    @property
    def vmean(self):
        return self.sum_volumes / self.nr_living

    @property
    def rcrit(self):
        return self.rmean / self.k

    @property
    def rhomin(self):
        return self.r[self.nr_living - 1] / self.r[0]

    @property
    def rhomax(self):
        return self.r[0] / self.r[0]

    @property
    def rhomean(self):
        return self.sum_radii / self.nr_living / self.r[0]

    @property
    def rhocrit(self):
        return -0.5 + np.sqrt(0.25 + self.rhomean / (self.k - self.rhomean))


    @property
    def thr(self):
        return float(self.nr_living) / float(self.nr_total)

    @property
    def vfrac(self):
        return self.sum_volumes / self.reference_volume

    @property
    def error_message(self):
        return capi.lib.rampr_get_error_message(self)
