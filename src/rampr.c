

#include "rampr.h"

//
// Distribution
// 

struct rampr_distribution_t * rampr_init_distribution ( const double * radii, unsigned nr_total ) {

    //
    // Memory allocation
    // 

    // Try to allocate the struct and check for success
    struct rampr_distribution_t * dist = malloc( sizeof(struct rampr_distribution_t) );
    if ( dist == NULL ) return NULL;

    // Try to allocate radius and volume array and check for success
    dist->radii   = malloc( sizeof(double) * nr_total );
    dist->volumes = malloc( sizeof(double) * nr_total );

    if ( dist->radii == NULL || dist->volumes == NULL ) {
        // The c standard 7.20.3.2/2 guarantees that calling free on
        // NULL is okay, so don't bother checking what exactly went wrong
        free( dist->radii ); free( dist->volumes ); free( dist );
        return NULL;
    }


    //
    // Setting member variables
    //

    // Copy radius array and sort it: dist->radii[0] is largest,
    // dist->radii[nr_total-1] is smallest radius
    for ( int i = 0; i < nr_total; i++ ) dist->radii[i] = radii[i];
    rampr_sort(dist->radii, nr_total);

    // Set total number of droplets, initialize counting of living/dead droplets
    dist->nr_total  = nr_total;
    dist->nr_living = nr_total;
    dist->nr_dead   = 0;

    // Count the number of non-evaporated droplets, set the radius of
    // evaporated ones to 0.
    for ( int i = nr_total - 1; i >= 0; i-- ) {
        if ( dist->radii[i] * dist->radii[i] * dist->radii[i] <  RAMPR_EVAPORATION_TOLERANCE ) {
            dist->radii[i] = 0.; dist->nr_living--; dist->nr_dead++;
        }
        else break; // Can do this because of the sorting above
    }

    // Initialize volumes, sum_radii, sum_volumes correctly
    dist->sum_radii   = 0.;
    dist->sum_volumes = 0.;

    for ( int i = 0; i < nr_total; i++ ) {
        dist->volumes[i] = dist->radii[i] * dist->radii[i] * dist->radii[i];

        dist->sum_radii   += dist->radii[i];
        dist->sum_volumes += dist->volumes[i];
    }

    return dist;
}


void rampr_free_distribution( struct rampr_distribution_t * dist ) {
    free( dist->radii );
    free( dist->volumes );
    free( dist );
}


struct rampr_distribution_t * rampr_merge_distributions ( const struct rampr_distribution_t * dist1,
                                                          const struct rampr_distribution_t * dist2,
                                                          bool remember_dead 
                                                        ) {
    double * radii;

    int nr1 = remember_dead ? dist1->nr_total : dist1->nr_living;
    int nr2 = remember_dead ? dist2->nr_total : dist2->nr_living;
    int nr  = nr1 + nr2;

    radii = malloc( sizeof(double) * nr );

    for ( int i = 0; i < nr1; i++ ) radii[i]     = dist1->radii[i];
    for ( int i = 0; i < nr2; i++ ) radii[nr1+i] = dist2->radii[i];

    return rampr_init_distribution(radii, nr);
}




//
// Simulation
//


struct rampr_simulation_t * rampr_init_simulation ( const double * radii, 
                                                    unsigned nr_total, 
                                                    double k_value 
                                                  ) {

    //
    // Memory allocation
    // 

    struct rampr_simulation_t * sim = malloc( sizeof(struct rampr_simulation_t) );
    if ( sim == NULL ) return NULL;

    //
    // Set member variables
    //

    sim->dist                = rampr_init_distribution(radii, nr_total);

    sim->time                = 0.;
    sim->k_value             = k_value;
    sim->volume_growth_rate  = (k_value - 1.) * 3 * sim->dist->nr_living;

    sim->error_tolerance     = RAMPR_DEFAULT_NORMED_ERROR_TOLERANCE;
    sim->nr_failed_adaptions = 0;
    sim->nr_steps            = 0;
    sim->target_dt           = RAMPR_DEFAULT_TARGET_DT;
    sim->last_dt             = 0;

    sim->reference_volume    = sim->dist->sum_volumes;

    // Set initial error message to empty string (null-char-terminated)
    sim->error_message[0]  = '\0';

    return sim;
}


void rampr_free_simulation ( struct rampr_simulation_t * sim ) {
    rampr_free_distribution( sim->dist );
    free( sim );
}




bool rampr_revive_droplets ( struct rampr_simulation_t * sim,
                             const double * radii,
                             unsigned n,
                             bool preserve_vfrac,
                             bool preserve_k_value
                           ) {

    //
    // Look how many given droplets (in radii) are alive
    //

    // Create a copy of radii that we can sort without side effects
    double * radii_ = malloc( sizeof(double) * n );
    if ( radii_ == NULL ) {
        sprintf( sim->error_message,
                 "Could not allocate temporary memory for reviving at t = %f",
                 sim->time
               );
        return false;
    }
    
    for ( int i = 0; i < n; i++ ) radii_[i] = radii[i];

    rampr_sort(radii_, n);

    unsigned nr_revive = n;

    // Count the number of living droplets in radii_
    for ( int i = n - 1; i >= 0; i-- ) {
        if ( radii_[i] * radii_[i] * radii_[i] <  RAMPR_EVAPORATION_TOLERANCE )
            nr_revive--;
        else
            break;
    }

    // Can't revive more droplets than have died (use 'rampr_add_droplets'
    // for that)
    if ( nr_revive > sim->dist->nr_dead ) {
        sprintf( sim->error_message,
                 "Tried to revive more particles (%d) than are dead (%d) at t = %f",
                 nr_revive, sim->dist->nr_dead, sim->time
               );
        free( radii_ );
        return false;
    }


    //
    // Revive dead droplets
    //
    
    double radius, volume;
    double sum_volumes_original = sim->dist->sum_volumes;
    int nr_living_original   = sim->dist->nr_living;

    for ( int i = 0; i < nr_revive; i++ ) {
        radius = radii_[i];
        volume = radii_[i] * radii_[i] * radii_[i];

        sim->dist->sum_radii    += radius;
        sim->dist->sum_volumes  += volume;


        sim->dist->radii[nr_living_original + i]   = radius;
        sim->dist->volumes[nr_living_original + i] = volume;

        sim->dist->nr_living++;
        sim->dist->nr_dead--;
    }

    //
    // Sort, update, and return
    //

    rampr_sort(sim->dist->radii,   sim->dist->nr_total);
    rampr_sort(sim->dist->volumes, sim->dist->nr_total);

    if ( preserve_vfrac )
        sim->reference_volume *= sim->dist->sum_volumes / sum_volumes_original;
    if ( preserve_k_value ) 
        sim->volume_growth_rate *= ( (double) sim->dist->nr_living ) / 
                                   ( sim->dist->nr_living - nr_revive );

    sim->k_value = 1 + sim->volume_growth_rate / (3. * sim->dist->nr_living);



    free( radii_ );

    return true;
}



bool rampr_add_droplets( struct rampr_simulation_t * sim, 
                         const double * radii, 
                         unsigned n, 
                         bool preserve_vfrac, 
                         bool preserve_k_value, 
                         bool include_dead 
                       ) {

    //
    // Combine sim->radii and radii, use this to get new distribution
    //

    double nr_total = n + sim->dist->nr_total;

    double * radii_ = malloc( sizeof(double) * ( nr_total ) );
    if ( radii_ == NULL ) {
        sprintf( sim->error_message,
                 "Could not allocate temporary memory for a new distribution at t = %f",
                 sim->time
               );
        return false;
    }


    for ( int i = 0; i < n; i++)
        radii_[i] = radii[i];
    for ( int i = 0; i < sim->dist->nr_total; i++)
        radii_[n + i] = sim->dist->radii[i];

    rampr_sort(radii_, nr_total);


    //
    // Change to new distribution
    //
    
    // Adapt nr_total if dead particles shall be left out
    if ( !include_dead ) {
        for ( int i = nr_total - 1; i >= 0; i-- ) {
            if ( radii_[i] * radii_[i] * radii_[i] < RAMPR_EVAPORATION_TOLERANCE ) 
                nr_total--;
            else 
                break;
        }

    }
    
    // Temporary backup of old distribution
    struct rampr_distribution_t * dist_before = sim->dist;

    // Create new one
    sim->dist = rampr_init_distribution(radii_, nr_total);
    if ( sim->dist == NULL ) {
        sprintf( sim->error_message,
                 "New distribution could not be allocated while adding droplets."
               );
        // Reattach the old distribution, such that calls that don't check
        // the status of the simulation won't result in segfaults
        sim->dist = dist_before;
        free( radii_ );
        return false;
    }


    //
    // Preserve the k_value and/or vfrac if desired
    //
    
    if ( preserve_vfrac )
        sim->reference_volume   *= sim->dist->sum_volumes / dist_before->sum_volumes;
    if ( preserve_k_value )
        sim->volume_growth_rate *= ( (double) sim->dist->nr_living ) / dist_before->nr_living;

    sim->k_value = 1 + sim->volume_growth_rate / (3. * sim->dist->nr_living);

    free( dist_before );

    return true;
}



//
// Helper functions
//

double rampr_third_root( const double v ) {
    return cbrt(v);
}


int rampr_cmpfkt( const void * a, const void * b ) {
    double dif = *((double *) b) - *((double*) a);
    return (dif > 0) - (dif < 0);
}

void rampr_sort( double * radii, unsigned nr_total ) {
    qsort(radii, nr_total, sizeof(double), rampr_cmpfkt);
}



double rampr_linear_time_until_evaporation( const struct rampr_simulation_t * sim, 
                                            double limit
                                          ) {

    double r_mean = sim->dist->sum_radii / sim->dist->nr_living; 
    double smallest_radius = sim->dist->radii[sim->dist->nr_living - 1];
    double smallest_volume = sim->dist->volumes[sim->dist->nr_living - 1];

    // this value should be ~ -3 if we look at particles evaporating soon
    double dvdt = 3 * ( sim->k_value * smallest_radius / r_mean - 1 );
    // the time derivative must be negative, if the droplet is to evaporate
    // at all (in linear approximation)
    if ( dvdt >= 0 ) return limit;
    // get the remaining lifetime of the smallest droplet  
    double delta_time = - smallest_volume / dvdt;
    // but return at most the limit value
    return delta_time < limit ? delta_time : limit;
}

double rampr_min_time_until_evaporation( const struct rampr_simulation_t * sim, 
                                         double limit
                                       ) {

    double r_mean = sim->dist->sum_radii / sim->dist->nr_living;
    double smallest_radius = sim->dist->radii[sim->dist->nr_living - 1];
    double smallest_volume = sim->dist->volumes[sim->dist->nr_living - 1];

    // this value should be ~ -3 if we look at particles evaporating soon
    double dvdt = 3 * ( sim->k_value * smallest_radius / r_mean - 1 );
    /*// the time derivative must be negative, if the droplet is to evaporate*/
    /*// at all (in linear approximation)*/
    if ( dvdt >= 0 ) return limit;
    // get the remaining lifetime of the smallest droplet  
    double delta_time = smallest_volume / 3.;
    // but return at most the limit value
    return delta_time < limit ? delta_time : limit;
}


double rampr_min_forward_time_until_evaporation( const struct rampr_simulation_t * sim,
                                                 double limit
                                               ) {

    double r_mean = sim->dist->sum_radii / sim->dist->nr_living;
    double smallest_radius = sim->dist->radii[sim->dist->nr_living - 1];

    // some helper variables
    double a = sim->k_value / r_mean;
    double b = a * smallest_radius;

    // check if the particle is growing (then return limit) or shrinking
    double dvdt = 3 * ( b - 1 );
    if ( dvdt >= 0 ) return limit;

    // get the remaining lifetime of the smallest droplet  
    /*double delta_time = - (log(1 - b) + b + b*b/2.) / (a*a*a);*/
    double delta_time = 0;
    for ( int k = 3; k <= 10; k++ ) {
        double power = 1.;
        for ( int i = 0; i < k; i++ )
            power *= b;
        delta_time += power / k;
    }

    delta_time /= a*a*a;

    // return at most the limit value
    return delta_time < limit ? delta_time : limit;
}


// check simulations for validity: has an error occured?
// this is checked before every integration step
bool rampr_is_valid( const struct rampr_simulation_t * sim ) {
    return sim->error_message[0] == '\0' ? true : false;
}

