

#include "rampr.h"



//
// Euler integration routine to calculate dvout
//

void rampr_euler_routine( double dt, const double * vin, const double * rin, 
                          double rin_mean, double k, unsigned nr_living,
                          double * dvout, double * verr, 
                          double * rbuffer, double ** vbuffers 
                        ) {

    //
    // Calculate the derivatives
    //

    rampr_volume_time_derivative( rin, rin_mean, k, nr_living, dvout );

    //
    // Euler method: dV is dV/dt * dt
    //

    for ( int i = 0; i < nr_living; i++ ) {
        dvout[i] = dvout[i] * dt;
    }

    // no error estimation here, so done
}


//
// Embedded Runge-Kutta method with Cash-Karp coefficients
// suggested in chapter 16.2 of 'numerical recipes in c'
//

void rampr_rkck_routine( double dt, const double * vin, const double * rin,
                         double rin_mean, double k, unsigned nr_living,
                         double * dvout, double * verr,
                         double * rbuffer, double ** vbuffers
                       ) {

    //
    // Define the Cash-Karp coefficients
    // Use static keyword so that they are not defined for every function call
    //

    static const double b21 = 0.2          , 
                        b31 = 3./40.       , b32 = 9./40.     ,  
                        b41 = 0.3          , b42 = -0.9       ,  b43 = 1.2         ,  
                        b51 = -11./54.     , b52 = 2.5        ,  b53 = -70./27.    , b54 = 35./27.        , 
                        b61 = 1631./55296. , b62 = 175./512.  ,  b63 = 575./13824. , b64 = 44275./110592. , b65 = 253./4096. ; 

    static const double c1  = 37./378.     , c3  = 250./621.  ,  c4  = 125./594.   , c6  = 512./1771. ; 

    static const double dc1 = c1 - 2825./27648.  , dc3 = c3 - 18575./48384.,
                        dc4 = c4 - 13525./55296. , dc5 =    - 277./14336., 
                        dc6 = c6 - 0.25          ;


    double * ak1      = vbuffers[0];
    double * ak2      = vbuffers[1];
    double * ak3      = vbuffers[2];
    double * ak4      = vbuffers[3];
    double * ak5      = vbuffers[4];
    double * ak6      = vbuffers[5];
    double * rtmp     = rbuffer;
    double   rmeantmp = rin_mean;


    //
    // 1. first step
    //

    rampr_volume_time_derivative( rin, rmeantmp, k, nr_living, ak1 );
    rmeantmp = 0;
    for ( int i = 0; i < nr_living; i++) {
        rtmp[i]   = rampr_third_root(vin[i] + b21 * ak1[i] * dt);
        rmeantmp += rtmp[i];
    }
    rmeantmp /= nr_living;


    //
    // 2. second step
    //

    rampr_volume_time_derivative( rtmp, rmeantmp, k, nr_living, ak2 );
    rmeantmp = 0;
    for ( int i = 0; i < nr_living; i++) {
        rtmp[i]   = rampr_third_root(vin[i] + ( b31 * ak1[i] + b32 * ak2[i] ) * dt);
        rmeantmp += rtmp[i];
    }
    rmeantmp /= nr_living;


    //
    // 3. third step
    //

    rampr_volume_time_derivative( rtmp, rmeantmp, k, nr_living, ak3 );
    rmeantmp = 0;
    for ( int i = 0; i < nr_living; i++) {
        rtmp[i]   = rampr_third_root(vin[i] + ( b41 * ak1[i] + b42 * ak2[i] + b43 * ak3[i] ) * dt);
        rmeantmp += rtmp[i];
    }
    rmeantmp /= nr_living;


    //
    // 4. fourth step
    //

    rampr_volume_time_derivative( rtmp, rmeantmp, k, nr_living, ak4 );
    rmeantmp = 0;
    for ( int i = 0; i < nr_living; i++) {
        rtmp[i]   = rampr_third_root(vin[i] + ( b51 * ak1[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i] ) * dt);
        rmeantmp += rtmp[i];
    }
    rmeantmp /= nr_living;

    //
    // 5. fifth step
    //

    rampr_volume_time_derivative( rtmp, rmeantmp, k, nr_living, ak5 );
    rmeantmp = 0;
    for ( int i = 0; i < nr_living; i++) {
        rtmp[i]   = rampr_third_root(vin[i] + ( b61 * ak1[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i] ) * dt);
        rmeantmp += rtmp[i];
    }
    rmeantmp /= nr_living;

    //
    // 6. sixth step
    //

    rampr_volume_time_derivative( rtmp, rmeantmp, k, nr_living, ak6 );

    // Get final output differences
    for ( int i = 0; i < nr_living; i++ ) {
        dvout[i] = ( c1 * ak1[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i] ) * dt;
    }
    // Estimate error between the fourth and fifth order method
    for ( int i = 0; i < nr_living; i++ ) {
        verr[i] = fabs( dc1 * ak1[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i] ) * dt;
    }
    
}

