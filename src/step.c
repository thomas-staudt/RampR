

#include "rampr.h"


//
// Evaluation of the volume time derivatives dvdt for given values of r,
// rmean and k
//

void rampr_volume_time_derivative ( const double * r, double r_mean, 
                                    double k, unsigned nr_living, double * dvdt
                                  ) {

   for ( int i = 0; i < nr_living; i++ )
       dvdt[i] = 3 * ( k * r[i] / r_mean - 1 );
}



//
// Every simulation step ends by calling this function
// It makes sure particles with volumes below threshold count as
// evaporated and it updates simulation parameters for the next step
//

int rampr_finish_simulation_step(struct rampr_simulation_t * sim) {


    //
    // Iterate over droplets, starting with the smallest
    // Check if the droplet are evaporated or not
    //
    
    int nr_evaporated = 0;
    int i = sim->dist->nr_living - 1;

    while ( sim->dist->volumes[i] < RAMPR_EVAPORATION_TOLERANCE ) {
        sim->dist->volumes[i] = 0.;
        sim->dist->radii[i]   = 0.;
        sim->dist->nr_living    --;
        sim->dist->nr_dead      ++;
        nr_evaporated           ++;
        i                       --;
    }

    //
    // Check if new droplet number is still high enough
    //

    if ( sim->dist->nr_living < RAMPR_MINIMAL_LIVING_DROPLETS ) {
        sprintf( sim->error_message, 
                 "Only %d droplets alive. Allowed minimum "
                 "number of living droplets is %d", 
                 sim->dist->nr_living, RAMPR_MINIMAL_LIVING_DROPLETS
               );
        return -1;
    }


    //
    // Update necessary variables and return
    //

    // k value
    sim->k_value = 1 + sim->volume_growth_rate / (3 * sim->dist->nr_living);

    // radius values
    sim->dist->sum_radii = 0;
    for ( i = 0; i < sim->dist->nr_living; i++ ) {
        sim->dist->radii[i]   = rampr_third_root(sim->dist->volumes[i]);
        sim->dist->sum_radii += sim->dist->radii[i];
    }

    sim->nr_steps++;
    
    return nr_evaporated;
}



//
// Make an euler step without dynamic adaption
// dt is exponentially increasing
//

int rampr_euler_step( struct rampr_simulation_t * sim, double max_dt, 
                      double * dvout, double * verr, double * rbuffer, double ** vbuffers
                    ) {

    //
    // dt mustn't be larger than time to next evaporation and max_dt
    //
    
    double dt = rampr_linear_time_until_evaporation(sim, sim->target_dt);
    if ( dt > max_dt ) dt = max_dt;

    //
    // Call the routine
    //

    double r_mean = sim->dist->sum_radii / sim->dist->nr_living; 

    rampr_euler_routine( dt, sim->dist->volumes, sim->dist->radii, r_mean, 
                         sim->k_value, sim->dist->nr_living, dvout, verr, rbuffer, vbuffers );


    //
    // Now, normally: check if to accept or reject the routine result.
    // For this Euler step we always accept, but give warning if
    // total volume expected and obtained by simulation differ to much
    //

    double total_volume_incr = 0;
    double target_volume_incr = dt * sim->volume_growth_rate;

    for ( int i = 0; i < sim->dist->nr_living; i++ ) total_volume_incr += dvout[i];

    if ( fabs(total_volume_incr - target_volume_incr) / sim->dist->sum_volumes > 
         RAMPR_INCR_TOLERANCE * (RAMPR_VOLUME_INCR_OFFSET + fabs(target_volume_incr)) * 
         (dt + RAMPR_DT_TOLERANCE_OFFSET)

       ) {
        fprintf( stderr, "Warning: Targeted (%.15f) and actual (%.15f) volume "
                         "increments differ strongly for t = %.15f\n", 
                 target_volume_incr, total_volume_incr, sim->time
               );
    }

    //
    // Update variables and finish the step
    //

    for ( int i = 0; i < sim->dist->nr_living; i++ ) {
        sim->dist->volumes[i]   += dvout[i];
    }

    sim->time              += dt;
    sim->last_dt            = dt;
    sim->target_dt         *= RAMPR_DEFAULT_VOLUME_INCR_FACTOR;
    sim->dist->sum_volumes += total_volume_incr;
    
    return rampr_finish_simulation_step(sim);
}


//
// Make an rkck step with dynamic timestep adaption
//

int rampr_rkck_step( struct rampr_simulation_t * sim, double max_dt, 
                     double * dvout, double * verr, double * rbuffer, double ** vbuffers
                   ) {
    
    // 
    // make sure that no negative volumes occur in rkck routine by making
    // dt small enough
    //

    /*double dt = rampr_min_forward_time_until_evaporation(sim, sim->target_dt);*/
    double dt = rampr_linear_time_until_evaporation(sim, sim->target_dt);
    if ( dt > max_dt ) dt = max_dt;

    //
    // apply the rkck routine until step is accepted (adaptive dt)
    //

    unsigned tries = 1; // count the tries
    double r_mean  = sim->dist->sum_radii / sim->dist->nr_living; 

    while ( true ) {

        //
        // make the calculation step and obtain the error estimates
        //

        rampr_rkck_routine( dt, sim->dist->volumes, sim->dist->radii, r_mean, 
                            sim->k_value, sim->dist->nr_living, dvout, verr, rbuffer, vbuffers );

        double normed_error = 0;
        for ( int i = 0; i < sim->dist->nr_living; i++ ) {
            if ( verr[i] > normed_error )
                normed_error = verr[i];
        }

        normed_error /= dt; // maximal error per unit time


        //
        // now decide if the step should be repeated, or if it should be
        // accepted; in both cases estimates for the new dt are needed
        // 
        
        //
        // case 1: repeat step with lower dt
        //
        
        if ( normed_error > sim->error_tolerance ) {

            // give up if too many tries are needed
            if ( tries > RAMPR_MAXIMAL_ADAPTION_TRIES ) {
                sprintf( sim->error_message, 
                         "Could not acquire sufficient accuracy (%e) with "
                         "maximal number of adaption steps (%d) at t = %.10f.",
                         sim->error_tolerance, RAMPR_MAXIMAL_ADAPTION_TRIES, sim->time
                       );

                return -1;
            }

            // error too high, so adapt dt
            dt = dt * 1.0 * pow(sim->error_tolerance / normed_error, 0.25);

            // check if the reduced dt is still reasonable
            if ( dt < RAMPR_MINIMAL_ACCEPTED_DT ) {
                sprintf( sim->error_message, 
                         "Adapted dt (%e) fell below minimal accepted dt (%e) at t = %.10f.",
                         dt, RAMPR_MINIMAL_ACCEPTED_DT, sim->time
                       );

                return -1;
            }

            // increment counters and repeat the while loop
            sim->nr_failed_adaptions++;
            tries++;
        }

        // 
        // case 2: estimate new dt for accepted step
        //

        else {

            // set new target_dt if 
            //     - well below tolerance
            if ( 5. * normed_error < sim->error_tolerance ) {
                double incr = pow(sim->error_tolerance / normed_error / 5., 0.2);
                sim->target_dt = dt * fmin(incr, RAMPR_MAXIMAL_DT_INCREMENT_FACTOR);
            }
            else
                sim->target_dt = dt;
            
            // break the loop and continue
            break;
        }
    }
    
    //
    // consistency check: targeted volume increment == simulated volume increment?
    //

    double target_volume_incr = dt * sim->volume_growth_rate;
    double total_volume_incr  = 0;

    for ( int i = 0; i < sim->dist->nr_living; i++ ) 
        total_volume_incr += dvout[i];

    if ( fabs(total_volume_incr - target_volume_incr) / sim->dist->sum_volumes > 
         RAMPR_INCR_TOLERANCE * (RAMPR_VOLUME_INCR_OFFSET + target_volume_incr) *
         (dt + RAMPR_DT_TOLERANCE_OFFSET)
       ) {
        fprintf( stderr, "Warning: Targeted (%.15f) and actual (%.15f) volume "
                         "increments differ for t = %.15f\n", 
                 target_volume_incr, total_volume_incr, sim->time
               );
    }

    //
    // update volume values and other parameters
    //

    for ( int i = 0; i < sim->dist->nr_living; i++ ) {
        sim->dist->volumes[i]  += dvout[i];
    }

    sim->time              += dt;
    sim->last_dt            = dt;
    sim->dist->sum_volumes += total_volume_incr;
    
    // 
    // check for evaporated particles, update k, xi, and r
    // 
    
    return rampr_finish_simulation_step(sim);
}


//
// Make an rkck step without dynamic adaption
// dt is exponentially increasing
//

int rampr_rkck_step_unadapted( struct rampr_simulation_t * sim, double max_dt, 
                               double * dvout, double * verr, double * rbuffer, double ** vbuffers
                             ) {
    
    //
    // dt mustn't be larger than time to next evaporation
    //

    double dt = rampr_min_forward_time_until_evaporation(sim, sim->target_dt);
    if ( dt > max_dt ) dt = max_dt;

    //
    // Call the routine
    //
    
    double r_mean = sim->dist->sum_radii / sim->dist->nr_living; 

    rampr_rkck_routine( dt, sim->dist->volumes, sim->dist->radii, r_mean, 
                        sim->k_value, sim->dist->nr_living, dvout, verr, rbuffer, vbuffers );


    //
    // Now, normally: check if to accept or reject the routine result.
    // For this Euler step we always accept, but give warning if
    // total volume expected and obtained by simulation differ to much
    //

    double target_volume_incr = dt * sim->volume_growth_rate;
    double total_volume_incr = 0;

    for ( int i = 0; i < sim->dist->nr_living; i++ ) total_volume_incr += dvout[i];

    if ( fabs(total_volume_incr - target_volume_incr) / sim->dist->sum_volumes > 
         RAMPR_INCR_TOLERANCE * (RAMPR_VOLUME_INCR_OFFSET + target_volume_incr) *
         (dt + RAMPR_DT_TOLERANCE_OFFSET)
       ) {
        fprintf( stderr, "Warning: Targeted (%.15f) and actual (%.15f) volume "
                         "increments differ strongly for t = %.15f\n", 
                 target_volume_incr, total_volume_incr, sim->time
               );
    }

    //
    // Update variables and finish the step
    //

    for ( int i = 0; i < sim->dist->nr_living; i++ ) {
        sim->dist->volumes[i]   += dvout[i];
    }

    sim->time              += dt;
    sim->last_dt            = dt;
    sim->target_dt         *= RAMPR_DEFAULT_VOLUME_INCR_FACTOR;
    sim->dist->sum_volumes += total_volume_incr;
    
    return rampr_finish_simulation_step(sim);
}
