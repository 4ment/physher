/*
 *  boot.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 6/12/2013.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "boot.h"

#include <math.h>
#include <stdio.h>

#include "gaussian.h"
#include "descriptivestats.h"
#include "utils.h"
#include "matrix.h"
#include "statistics.h"

// Normal interval
double * boot_ci_norm( const double *thetas, int n, double theta_hat, double alpha ){
    double pp = 1.0-alpha;
    double qlb = (1.0-pp)/2.0;
    double qub = (1.0-pp)/2.0+pp;
    double mean = dmean(thetas, n);
    double s = standard_deviation(thetas, n, mean);
    double *ci = dvector(2);
    ci[0] = theta_hat + qnorm(qlb,0,1)*s;
    ci[1] = theta_hat + qnorm(qub/2.0,0,1)*s;
    return ci;
}

// Normal interval
double boot_ci_norm2( const double *thetas, int n, double theta_hat, double q ){
    double mean = dmean(thetas, n);
    double s = standard_deviation(thetas, n, mean);
    return theta_hat + qnorm(q,0,1)*s;
}


//Bias correction
static double _bias_correction( const double *thetas, int b, double theta_hat ){
    double n = 0;
    for ( int i = 0; i < b; i++ ) {
        if( thetas[i] < theta_hat ){
            n++;
        }
    }
    n /= b;
    return qnorm(n, 0, 1);
}

// non-parametric estimation of the acceleration
static double _acceleration( double *thetai, int n, double theta_hat ){
    double a = 0;
    double sq = 0;
    double cb = 0;
    
    for ( int i = 0; i < n; i++ ) {
        double diff = theta_hat-thetai[i];
        sq += diff * diff;
        cb += diff * diff * diff;
    }
    a = cb/(6.0*pow(sq,1.5));
    return a;
}

static double _acceleration_weighted( double *thetai, int n, double theta_hat, const double *weights ){
    double a = 0;
    double sq = 0;
    double cb = 0;
    
    for ( int i = 0; i < n; i++ ) {
        double diff = theta_hat-thetai[i];
        sq += diff * diff * weights[i];
        cb += diff * diff * diff * weights[i];
    }
    a = cb/(6.0*pow(sq,1.5));
    return a;
}

// alpha is the endpoint (e.g. 0.025 or 0.975 for 95% CI)
static double _bias_corrected_endpoint( double z0, double a, double alpha ){
    return pnorm(z0 + (z0+qnorm(alpha,0,1))/(1.0-a*(z0+qnorm(alpha,0,1)) ),0,1);
}

// thetai is the jackknifed estimate of the parameter, obtained by omitting observation i from the data. Not necessarily ordered.
// n: length of thetai vector
// thetab: estimates of the bootstrap
// b: number of bootstrap
// theta_hat: point estimate of the parameter of interest obtained from the original dataset.
// q: quantile
double bootci_BCa(double *thetai, int n, const double *thetab, int b, double theta_hat, double q){
    double *v2 = clone_dvector( thetai, n );
    qsort(v2, n, sizeof(double), qsort_asc_dvector);
    if( n == 1 || v2[0] == v2[n-1]){
        double v = v2[0];
        free(v2);
        return v;
    }
    free(v2);
    
    double z0 = _bias_correction(thetab, b, theta_hat);
    double a = _acceleration(thetai, n, theta_hat);
    double qc = _bias_corrected_endpoint(z0, a, q);
    
    if( qc == 0 ){
        return thetab[0];
    }
    
    //fprintf(stderr, "z0: %e theta_hat: %e\n", z0, theta_hat);
    //fprintf(stderr, "a: %e\n", a);
    //fprintf(stderr, "qc: %e\n", qc);
    return dpercentile(thetab, b, qc);
}

double bootci_BCa_weighted(double *thetai, int n, const double *weights, const double *thetab, int b, double theta_hat, double q){
    double *v2 = clone_dvector( thetai, n );
    qsort(v2, n, sizeof(double), qsort_asc_dvector);
    if( n == 1 || v2[0] == v2[n-1]){
        free(v2);
        return thetai[0];
    }
    free(v2);
    
    double z0 = _bias_correction(thetab, b, theta_hat);
    double a = _acceleration_weighted(thetai, n, theta_hat, weights);
    double qc = _bias_corrected_endpoint(z0, a, q);
    if( qc == 0 ){
        return thetab[0];
    }
    //    if(qc == 0.0 ){
    //        fprintf(stderr, "z0: %e theta_hat: %e\n", z0, theta_hat);
    //        fprintf(stderr, "a: %e\n", a);
    //        fprintf(stderr, "qc: %e\n", qc);
    //    }
    return dpercentile(thetab, b, qc);
}

double bootci_BCa_weighted_debug(double *thetai, int n, const double *weights, const double *thetab, int b, double theta_hat, double q){
    double *v2 = clone_dvector( thetai, n );
    qsort(v2, n, sizeof(double), qsort_asc_dvector);
    if( n == 1 || v2[0] == v2[n-1]){
        free(v2);
        return thetai[0];
    }
    free(v2);
    
    double z0 = _bias_correction(thetab, b, theta_hat);
    double a = _acceleration_weighted(thetai, n, theta_hat, weights);
    double qc = _bias_corrected_endpoint(z0, a, q);

    fprintf(stderr, "z0: %e theta_hat: %e\n", z0, theta_hat);
    fprintf(stderr, "a: %e\n", a);
    fprintf(stderr, "qc: %e\n", qc);
    
    print_dvector(thetai, n);
    print_dvector(thetab, b);
    
    
    for (int i = 0; i < n; i++ ) {
        for ( int j = 0; j < weights[i]; j++ ) {
            fprintf(stderr, "%e,", thetai[i]);
        }
    }
    fprintf(stderr, "\n\n");
    
    
    
    for (int i = 0; i < b; i++ ) {
        fprintf(stderr, "%e,", thetab[i]);
    }
    fprintf(stderr, "\n");

    return dpercentile(thetab, b, qc);
}

