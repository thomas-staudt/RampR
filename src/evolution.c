

#include "rampr.h"


// 
// Make exponentially increasing euler evolution steps until 'sim->time = until_time'
// The euler routine does not estimate discretization errors, so setting 'adapted'
// is without effect
// 

int rampr_euler_evolution_until( struct rampr_simulation_t * sim,
                                 double until_time,
                                 unsigned until_nr,
                                 bool adapted
                               ) {
    //
    // Check sim for validity (has an uncaught error occured before?)
    //

    if ( !rampr_is_valid(sim) ) return -1;


    //
    // Create buffers for the routine
    //
    
    // Volume differences
    double * dvout = malloc( sizeof(double) * sim->dist->nr_total);
    // These are not necessary for simple euler routine
    double *  verr     = NULL;
    double ** vbuffers = NULL; 
    double *  rbuffer  = NULL;


    //
    // Make as many steps as needed to arrive at 'until_time',
    // count the evaporated particles
    //
    
    int nr_evaporated = 0;
    int total_nr_evaporated = 0;

    double max_dt;

    while ( sim->time < until_time && sim->dist->nr_living > until_nr ) {

        max_dt = until_time - sim->time + RAMPR_DT_EXCESS;

        nr_evaporated = rampr_euler_step(sim, max_dt, dvout, verr, rbuffer, vbuffers);
        // Check for errors (negative nr_evaporated) and propagate them
        if ( nr_evaporated < 0 ) {
            free ( dvout );
            return -1;
        }

        // Add to total number and continue
        total_nr_evaporated += nr_evaporated;
    }

    // Cleanup
    free ( dvout );

    return total_nr_evaporated;
}



// 
// Make a simple euler evolution step
// Probably has bad performance, since memory is allocated/deallocated on
// each call
// 

int rampr_euler_evolution_single( struct rampr_simulation_t * sim,
                                  double until_time,
                                  bool adapted
                                ) {

    //
    // Check sim for validity (has an uncaught error occured before?)
    //

    if ( !rampr_is_valid(sim) ) return -1;

    //
    // Create buffers for the routine
    //
    
    // Volume differences
    double * dvout = malloc( sizeof(double) * sim->dist->nr_total);
    // These are not necessary for simple euler routine
    double *  verr     = NULL;
    double ** vbuffers = NULL; 
    double *  rbuffer  = NULL;

    //
    // Make the Euler step and return
    //

    // Maximal allowed dt for the next step
    double max_dt; 

    if ( until_time > 0 )
        max_dt = fmax(0., until_time - sim->time) + RAMPR_DT_EXCESS;
    else
        max_dt = RAMPR_MAXIMAL_DT; 

    int nr_evaporated = rampr_euler_step(sim, max_dt, dvout, verr, rbuffer, vbuffers);

    // Cleanup
    free ( dvout );

    return nr_evaporated;
}




// 
// Make Runge-Kutta evolution steps, using the rkck-routine, until sim->time = until_time
// 'adapted' decides if timesteps dt adapted to error estimates (Cash-Karp)
// are used, or if they are chosen exponentially increasing (as with the
// euler method)
// 

int rampr_rkck_evolution_until( struct rampr_simulation_t * sim,
                                double until_time,
                                unsigned until_nr,
                                bool adapted
                              ) {
    //
    // Check for validity (have uncaught errors occured?)
    //

    if ( !rampr_is_valid(sim) ) return -1;


    //
    // Create buffers for the routine
    //
    
    // rkck routine needs buffers for volume differences, volume errors,
    // radii, and 6 further evaluations of the time derivative
    double * dvout       = malloc( sizeof(double) * sim->dist->nr_total );
    double * verr        = malloc( sizeof(double) * sim->dist->nr_total );
    double * rbuffer     = malloc( sizeof(double) * sim->dist->nr_total );

    double * vbuffers[6];
    for ( int i = 0; i < 6; i++ )
        vbuffers[i] = malloc( sizeof(double) * sim->dist->nr_total );


    //
    // Make as many steps as needed to arrive at 'until_time',
    // count the number of evaporated particles
    //
    
    int nr_evaporated = 0;
    int total_nr_evaporated = 0;

    double max_dt = 0.;

    while ( sim->time < until_time && sim->dist->nr_living > until_nr  ) {

        max_dt = until_time - sim->time;

        if ( adapted )
            nr_evaporated = rampr_rkck_step(sim, max_dt, dvout, verr, rbuffer, vbuffers);
        else
            nr_evaporated = rampr_rkck_step_unadapted(sim, max_dt, dvout, verr, rbuffer, vbuffers);

        // Check for errors and propagate them
        if ( nr_evaporated < 0 ) {
            free ( dvout ); free( verr ); free( rbuffer );
            for ( int i = 0; i < 6; i++ ) free( vbuffers[i] );
            return -1;
        }

        // Add to total number and continue
        total_nr_evaporated += nr_evaporated;
    }

    // Cleanup
    free ( dvout ); free( verr ); free( rbuffer );
    for ( int i = 0; i < 6; i++ ) free( vbuffers[i] );


    return total_nr_evaporated;
}


//
// Make a single rkck step with option of adapted or non-adapted dt
// Probably very inefficient, since much memory is allocated/deallocated
//

int rampr_rkck_evolution_single( struct rampr_simulation_t * sim, double until_time, bool adapted ) {

    //
    // Check for validity (have uncaught errors occured?)
    //

    if ( !rampr_is_valid(sim) ) return -1;
    
    //
    // Create buffers for the routine
    //
    
    double * dvout       = malloc( sizeof(double) * sim->dist->nr_total );
    double * verr        = malloc( sizeof(double) * sim->dist->nr_total );
    double * rbuffer     = malloc( sizeof(double) * sim->dist->nr_total );

    double * vbuffers[6];
    for ( int i = 0; i < 6; i++ )
        vbuffers[i] = malloc( sizeof(double) * sim->dist->nr_total );


    //
    // Make the rkck step and return
    //

    // Maximal allowed dt for the next step
    double max_dt; 

    if ( until_time > 0 )
        max_dt = fmax(0., until_time - sim->time) + RAMPR_DT_EXCESS;
    else
        max_dt = RAMPR_MAXIMAL_DT; 


    int nr_evaporated = 0;

    if ( adapted )
        nr_evaporated = rampr_rkck_step(sim, max_dt, dvout, verr, rbuffer, vbuffers);
    else
        nr_evaporated = rampr_rkck_step_unadapted(sim, max_dt, dvout, verr, rbuffer, vbuffers);

    // Cleanup
    free ( dvout ); free( verr ); free( rbuffer );
    for ( int i = 0; i < 6; i++ ) free ( vbuffers[i] );

    return nr_evaporated;
}


