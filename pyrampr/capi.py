
# Makes the functions in librampr.so available in python

from os import path
from numpy import ctypeslib, float64
import ctypes

#
# Define types and variables
#

# Introduce the ctypes types needed
_uint      = ctypes.c_uint
_int       = ctypes.c_int
_double    = ctypes.c_double
_char_p    = ctypes.c_char_p
_void_p    = ctypes.c_void_p
_void      = None
_float_vec = ctypeslib.ndpointer( dtype=float64, 
                                  ndim=1, 
                                  flags='C_CONTIGUOUS'
                                )


# Load the c library librampr.so
module_dir = path.dirname(path.abspath(__file__))
lib_dir    = path.join(module_dir, '..', 'lib')
lib_name   = "librampr.so"
lib        = ctypeslib.load_library(lib_name, lib_dir)



lib.EVAPORATION_TOLERANCE = 1e-20

#
# Specify signature and return types of the c-functions
#


# init_simulation
lib.rampr_init_simulation.argtypes          = [ _float_vec, _uint, _double ]
lib.rampr_init_simulation.restype           = _void_p

# free_simulation
lib.rampr_free_simulation.argtypes          = [ _void_p ]
lib.rampr_free_simulation.restype           = _void

# evolve -- single steps
lib.rampr_euler_evolution_single.argtypes   = [ _void_p, _double, _int ]
lib.rampr_euler_evolution_single.restype    = _int

lib.rampr_rkck_evolution_single.argtypes    = [ _void_p, _double, _int ]
lib.rampr_rkck_evolution_single.restype     = _int

# evolve -- until a certain time value
lib.rampr_euler_evolution_until.argtypes    = [ _void_p, _double, _uint, _int ]
lib.rampr_euler_evolution_until.restype     = _int

lib.rampr_rkck_evolution_until.argtypes     = [ _void_p, _double, _uint, _int ]
lib.rampr_rkck_evolution_until.restype      = _int


# change the distribution
lib.rampr_revive_droplets.argtypes          = [ _void_p, _float_vec, _uint, _int, _int ]
lib.rampr_revive_droplets.restype           = _int

lib.rampr_add_droplets.argtypes             = [ _void_p, _float_vec, _uint, _int, _int, _int ]
lib.rampr_add_droplets.restype              = _int



# radii
lib.rampr_get_radii.argtypes                = [ _void_p ]
lib.rampr_get_radii.restype                 = ctypes.POINTER(_double)

lib.rampr_get_sum_radii.argtypes            = [ _void_p ]
lib.rampr_get_sum_radii.restype             = _double

# volumes
lib.rampr_get_volumes.argtypes              = [ _void_p ]
lib.rampr_get_volumes.restype               = ctypes.POINTER(_double)

lib.rampr_get_sum_volumes.argtypes          = [ _void_p ]
lib.rampr_get_sum_volumes.restype           = _double

# nr_initial
lib.rampr_get_nr_total.argtypes           = [ _void_p ]
lib.rampr_get_nr_total.restype            = _uint

# nr_living
lib.rampr_get_nr_living.argtypes            = [ _void_p ]
lib.rampr_get_nr_living.restype             = _uint

# nr_dead
lib.rampr_get_nr_dead.argtypes              = [ _void_p ]
lib.rampr_get_nr_dead.restype               = _uint

# time
lib.rampr_get_time.argtypes                 = [ _void_p ]
lib.rampr_get_time.restype                  = _double

lib.rampr_set_time.argtypes                 = [ _void_p, _double ]
lib.rampr_set_time.restype                  = _void

lib.rampr_get_last_dt.argtypes              = [ _void_p ]
lib.rampr_get_last_dt.restype               = _double

lib.rampr_get_target_dt.argtypes            = [ _void_p ]
lib.rampr_get_target_dt.restype             = _double

lib.rampr_set_target_dt.argtypes            = [ _void_p, _double ]
lib.rampr_set_target_dt.restype             = _void

# error_tolerance
lib.rampr_get_error_tolerance.argtypes      = [ _void_p ]
lib.rampr_get_error_tolerance.restype       = _double

lib.rampr_set_error_tolerance.argtypes      = [ _void_p, _double ]
lib.rampr_set_error_tolerance.restype       = _void

# nr_failed_adaptions
lib.rampr_get_nr_failed_adaptions.argtypes  = [ _void_p ]
lib.rampr_get_nr_failed_adaptions.restype   = _uint

# nr_steps
lib.rampr_get_nr_steps.argtypes             = [ _void_p ]
lib.rampr_get_nr_steps.restype              = _uint

# k_value
lib.rampr_get_k_value.argtypes              = [ _void_p ]
lib.rampr_get_k_value.restype               = _double

lib.rampr_set_k_value.argtypes              = [ _void_p, _double ]
lib.rampr_set_k_value.restype               = _void

# volume_growth_rate
lib.rampr_get_volume_growth_rate.argtypes   = [ _void_p ]
lib.rampr_get_volume_growth_rate.restype    = _double

lib.rampr_set_volume_growth_rate.argtypes   = [ _void_p, _double ]
lib.rampr_set_volume_growth_rate.restype    = _void

# reference_volume
lib.rampr_get_reference_volume.argtypes     = [ _void_p ]
lib.rampr_get_reference_volume.restype      = _double

lib.rampr_set_reference_volume.argtypes     = [ _void_p, _double ]
lib.rampr_set_reference_volume.restype      = _void

# error message
lib.rampr_get_error_message.argtypes        = [ _void_p ]
lib.rampr_get_error_message.restype         = _char_p
