

#include "rampr.h"

//
// These getter/setter functions provide an api to the rampr_simulation_t
// struct that can be used by other programming languages
//



//
// Getter functions
//

const double * rampr_get_volumes( const struct rampr_simulation_t * sim ) {
    return sim->dist->volumes;
}

const double * rampr_get_radii( const struct rampr_simulation_t * sim ) {
    return sim->dist->radii;
}

unsigned rampr_get_nr_total( const struct rampr_simulation_t * sim ) {
    return sim->dist->nr_total;
}

unsigned rampr_get_nr_living( const struct rampr_simulation_t * sim ) {
    return sim->dist->nr_living;
}

unsigned rampr_get_nr_dead( const struct rampr_simulation_t * sim ) {
    return sim->dist->nr_dead;
}

double rampr_get_sum_volumes( const struct rampr_simulation_t * sim ) {
    return sim->dist->sum_volumes;
}

double rampr_get_sum_radii( const struct rampr_simulation_t * sim ) {
    return sim->dist->sum_radii;
}

double rampr_get_time( const struct rampr_simulation_t * sim ) {
    return sim->time;
}

double rampr_get_last_dt( const struct rampr_simulation_t * sim ) {
    return sim->last_dt;
}

double rampr_get_error_tolerance( const struct rampr_simulation_t * sim ) {
    return sim->error_tolerance;
}

unsigned rampr_get_nr_failed_adaptions( const struct rampr_simulation_t * sim ) {
    return sim->nr_failed_adaptions;
}

unsigned rampr_get_nr_steps( const struct rampr_simulation_t * sim ) {
    return sim->nr_steps;
}

double rampr_get_target_dt( const struct rampr_simulation_t * sim ) {
    return sim->target_dt;
}

double rampr_get_k_value( const struct rampr_simulation_t * sim ) {
    return sim->k_value;
}

double rampr_get_volume_growth_rate( const struct rampr_simulation_t * sim ) {
    return sim->volume_growth_rate;
}

double rampr_get_reference_volume( const struct rampr_simulation_t * sim ) {
    return sim->reference_volume;
}


const char * rampr_get_error_message( const struct rampr_simulation_t * sim ) {
    return sim->error_message;
}




//
// Setter functions
//


void rampr_set_time( struct rampr_simulation_t * sim, double time ) {
    sim->time = time;
}

void rampr_set_k_value( struct rampr_simulation_t * sim, double k_value ) {
    sim->k_value = k_value;
    sim->volume_growth_rate = (k_value - 1) * 3 * sim->dist->nr_living;
}

void rampr_set_volume_growth_rate( struct rampr_simulation_t * sim, 
                                   double volume_growth_rate 
                                 ) {
    sim->volume_growth_rate = volume_growth_rate;
    sim->k_value = 1. + volume_growth_rate / (3. * sim->dist->nr_living);
}

void rampr_set_target_dt( struct rampr_simulation_t * sim, 
                          double target_dt 
                        ) {
    sim->target_dt = target_dt;
}


void rampr_set_error_tolerance( struct rampr_simulation_t * sim, 
                                double error_tolerance 
                              ) {
    sim->error_tolerance = error_tolerance;
}


void rampr_set_reference_volume( struct rampr_simulation_t * sim, 
                                 double reference_volume 
                              ) {
    sim->reference_volume = reference_volume;
}

