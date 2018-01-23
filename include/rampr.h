
//
// Header file for the rampr simulation library
//

#ifndef __RAMPR2
#define __RAMPR2


#include<math.h>
#include<stdio.h>
#include<stdlib.h>

//
// Constants
//

#define RAMPR_INCR_TOLERANCE                   1e-5
#define RAMPR_VOLUME_INCR_OFFSET               1e-5
#define RAMPR_DT_TOLERANCE_OFFSET              1e-10
#define RAMPR_EVAPORATION_TOLERANCE            1e-20

#define RAMPR_MAXIMAL_DT                       1e10
#define RAMPR_DT_EXCESS                        1e-10

#define RAMPR_MINIMAL_ACCEPTED_DT              1e-40
#define RAMPR_MAXIMAL_DT_INCREMENT_FACTOR      1e20

#define RAMPR_MINIMAL_LIVING_DROPLETS          10
#define RAMPR_MAXIMAL_ADAPTION_TRIES           50

#define RAMPR_DEFAULT_TARGET_DT                1e-6
#define RAMPR_DEFAULT_NORMED_ERROR_TOLERANCE   1e-6 
#define RAMPR_DEFAULT_VOLUME_INCR_FACTOR       1.000002302587744 
                                               // = 10^(1e-6)



// Some syntactic sugar
typedef int bool;

#define true 1
#define false 0



// 
// The distribution struct
//


// A collection of droplets
struct rampr_distribution_t {

    // Radius and volume values for the droplets
    // are ordered such that they are shrinking in size
    // with growing index, i.e. radii[0] >= radii[1] >= ...
    // Radii and volumes are connected via radius^3 = volume
    double * volumes;
    double * radii;

    // Integrals over volumes and radii
    double sum_volumes;
    double sum_radii;


    // Total number of particles in the distribution
    // is the initial number of droplets by default
    unsigned nr_total;
    // Number of particles which are non-zero
    unsigned nr_living;
    // Number of particles which are evaporated
    unsigned nr_dead;

};

//
// The simulation struct
//
// Contains all necessary information to make an integration step on the
// volume values
// Running a simulation consists of the following steps:
//
//      * Initialize a simulation by providing a radius/volume array and
//        a initial k/xi value.
//      * Use a 'evolve'-function on the simulation, which will evolve the
//        distribution and simulation parameters in time
//

struct rampr_simulation_t {

    // Current distribution
    struct rampr_distribution_t * dist;

    // Memory for adaptive step size drivers
    double   target_dt;
    double   last_dt;
    double   error_tolerance;
    unsigned nr_failed_adaptions;
    unsigned nr_steps;

    // Simulation parameters
    double time;
    double k_value;
    double volume_growth_rate;

    // Reference value for sum_volumes, normally the 
    // initial value, except particles are manually 
    // "revived" or modified in another way 
    double reference_volume;

    // Somewhere to place error messages
    char error_message[1024];

};



//
// Basic methods for rampr_distribution_t
//

// Initialization and destruction of the struct
struct rampr_distribution_t * rampr_init_distribution   ( const double * radii, unsigned nr_initial );
void                          rampr_free_distribution   ( struct rampr_distribution_t * dist );

// Take two distributions and create a new one that contains the droplets
// of both dist1 and dist2, with an option to remember or discard dead
// droplets
struct rampr_distribution_t * rampr_merge_distributions  ( const struct rampr_distribution_t * dist1,
                                                           const struct rampr_distribution_t * dist2,
                                                           int remember_dead 
                                                         );


//bool                          rampr_refill_distribution ( struct rampr_distribution_t * dist,
                                                          //const double * fill_radii,
                                                          //unsigned n
                                                        //);


//
// Basic methods for rampr_simulation_t
//

struct rampr_simulation_t * rampr_init_simulation        ( const double * radii, unsigned nr_initial, double k_value );
void                        rampr_free_simulation        ( struct rampr_simulation_t * sim );



void                        rampr_volume_time_derivative ( const double * r, double r_mean, double k, unsigned nr_living, double * dvdt );



void                        rampr_euler_routine          ( double dt, const double * vin, const double * rin, double rin_mean, double k, 
                                                           unsigned nr_living, double * vout, double * verr, double * rbuffer, double ** vbuffers );
void                        rampr_rkck_routine           ( double dt, const double * vin, const double * rin, double rin_mean, double k, 
                                                           unsigned nr_living, double * vout, double * verr, double * rbuffer, double ** vbuffers );


int                         rampr_euler_step             ( struct rampr_simulation_t * sim, double max_dt, 
                                                           double * dvout, double * verr, double * rbuffer, double ** vbuffers );
int                         rampr_rkck_step              ( struct rampr_simulation_t * sim, double max_dt,
                                                           double * dvout, double * verr, double * rbuffer, double ** vbuffers );
int                         rampr_rkck_step_unadapted    ( struct rampr_simulation_t * sim, double max_dt,
                                                           double * dvout, double * verr, double * rbuffer, double ** vbuffers );


int                         rampr_euler_evolution_until  ( struct rampr_simulation_t * sim, double until_time, 
                                                           unsigned until_nr, bool adapted );

int                         rampr_rkck_evolution_until   ( struct rampr_simulation_t * sim, double until_time, 
                                                           unsigned until_nr, bool adapted );

int                         rampr_euler_evolution_single ( struct rampr_simulation_t * sim, double until_time, bool adapted );

int                         rampr_rkck_evolution_single  ( struct rampr_simulation_t * sim, double until_time, bool adapted );



bool                        rampr_revive_droplets        ( struct rampr_simulation_t * sim, const double * radii, unsigned n, 
                                                           bool preserve_vfrac, bool preserve_k_value );

bool                        rampr_add_droplets           ( struct rampr_simulation_t * dist, const double * radii, unsigned n, 
                                                           bool preserve_vfrac, bool preserve_k_value, 
                                                           bool include_dead );

//
// Setter and getter functions for simulations, needed for the pyrampr wrapper
//

// Get
const double *              rampr_get_volumes            ( const struct rampr_simulation_t * sim );
const double *              rampr_get_radii              ( const struct rampr_simulation_t * sim );
unsigned                    rampr_get_nr_initial         ( const struct rampr_simulation_t * sim );
unsigned                    rampr_get_nr_living          ( const struct rampr_simulation_t * sim );
unsigned                    rampr_get_nr_dead            ( const struct rampr_simulation_t * sim );
double                      rampr_get_time               ( const struct rampr_simulation_t * sim );
double                      rampr_get_last_dt            ( const struct rampr_simulation_t * sim );
double                      rampr_get_error_tolerance    ( const struct rampr_simulation_t * sim );
unsigned                    rampr_get_nr_failed_adaptions( const struct rampr_simulation_t * sim );
unsigned                    rampr_get_nr_steps           ( const struct rampr_simulation_t * sim );
double                      rampr_get_target_dt          ( const struct rampr_simulation_t * sim );
double                      rampr_get_sum_volumes        ( const struct rampr_simulation_t * sim );
double                      rampr_get_sum_radii          ( const struct rampr_simulation_t * sim );
double                      rampr_get_k_value            ( const struct rampr_simulation_t * sim );
double                      rampr_get_volume_growth_rate ( const struct rampr_simulation_t * sim );
double                      rampr_get_reference_volume   ( const struct rampr_simulation_t * sim );
const char *                rampr_get_error_message      ( const struct rampr_simulation_t * sim );

// Set
void                        rampr_set_time               ( struct rampr_simulation_t * sim, double time );
void                        rampr_set_k_value            ( struct rampr_simulation_t * sim, double k_value );
void                        rampr_set_volume_growth_rate ( struct rampr_simulation_t * sim, double volume_growth_rate );
void                        rampr_set_target_dt          ( struct rampr_simulation_t * sim, double target_dt );
void                        rampr_set_error_tolerance    ( struct rampr_simulation_t * sim, double error_tolerance );
void                        rampr_set_reference_volume   ( struct rampr_simulation_t * sim, double reference_volume );




//
// Other helper functions
//

double                      rampr_third_root             ( const double v );
int                         rampr_cmpfkt                 ( const void * a, const void * b );
void                        rampr_sort                   ( double * radii, unsigned nr_initial );


//
// Different ways to approximate the time until evaporation of the smallest
// droplet
// Return limit if this time is greater than the limit value
//

double                      rampr_linear_time_until_evaporation      ( const struct rampr_simulation_t * sim, double limit );
double                      rampr_min_time_until_evaporation         ( const struct rampr_simulation_t * sim, double limit );
double                      rampr_min_forward_time_until_evaporation ( const struct rampr_simulation_t * sim, double limit );


int                         rampr_finish_simulation_step ( struct rampr_simulation_t * sim );
bool                        rampr_is_valid               ( const struct rampr_simulation_t * sim );


#endif // __RAMPR2
