/*
 *  gamma.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/24/11.
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


#include "gamma.h"

#include <math.h>
#include <float.h>

#include "utils.h"
#include "gaussian.h"
#include "random.h"

#define ITMAX 100
#define FPMIN 1.0e-30

//#define EPS 3.0e-7

static double gser( double *gamser, double a, double x );
static double gcf( double *gammcf, double a, double x );


#pragma mark -
#pragma mark Gamma function

double dgamma(const double x, const double alpha, const double beta){
    return pow(beta, alpha)*pow(x,alpha-1.0)*exp(-beta*x)/gamm(alpha);
}

double dloggamma(const double x, const double alpha, const double beta){
    return log(pow(beta, alpha)*pow(x,alpha-1.0)) -beta*x - gammln(alpha);
}

// Cumulative distribution function of Gamma function
double pgamma( const double x, const double alpha, const double beta ){
	if (x < 0.) error("bad x in pgamma");
	return gammp(alpha,beta*x);
}

// Inverse cumulative distribution function of Gamma function
double qgamma( const double p, const double alpha, const double beta ){
	if (p < 0. || p >= 1.) error("bad p in pgamma");
	return invgammp(p,alpha)/beta;
}

double rgamma(double aa) {
	double a = aa;
	if(aa < 1.0){
		a += 1.0;
	}
	double d, c, x, v, u;
	d = a - 1. / 3.;
	c = 1. / sqrt(9. * d);
	for (;;) {
		do {
			x = rnorm();
			v = 1. + c * x;
		} while (v <= 0.);
		v = v * v * v;
		u = random_double();
		if (u < 1. - 0.0331 * (x * x) * (x * x)) {
			break;
		}
		if (log(u) < 0.5 * x * x + d * (1. - v + log(v))) {
			break;
		}
	}
	if(aa < 1.0){
		double u;
		do{
			u = random_double();
		} while(u == 0);
		return (d * v) * pow(u, 1.0/aa);
	}
	return (d * v);
}

#pragma mark -
#pragma mark Incomplete Gamma function

// Returns the incomplete gamma function P(a,x).
double gammp( const double a, const double x ){
	double gamser = 0;
	double gammcf;
	if (x < 0.0 || a < 0.0){
		fprintf(stderr, "Invalid arguments in routine gammp a=%f x=%f",a,x);
		exit(1);
	}
	if (x < (a+1.0)) {          // Use the series representation
		gser( &gamser, a, x );
		return gamser;
	}
	else {                      // Use the continued fraction representation
		gcf( &gammcf, a, x ); 
		return 1.0-gammcf;      //and take its complement.
	}
}


// Returns the incomplete gamma function P (a, x) evaluated by its series representation as gamser.
// Also returns ln Γ(a) as gln.
double gser( double *gamser, double a, double x){
	int n;
	double sum,del,ap;
	double gln = gammln(a);
	
	if (x <= 0.0) {
		if (x < 0.0) return NAN; // error("x less than 0 in routine gser");
		*gamser=0.0; 
		return gln;
	} 
	else {
		ap = a;
		del = sum = 1.0/a;
		for ( n = 0; n < ITMAX; n++ ) {
			++ap; 
			del *= x/ap;
			sum += del; 
			if (fabs(del) < fabs(sum)*DBL_EPSILON) {
				*gamser = sum * exp(-x+a*log(x)-gln); 
				return gln;
			}
		} 
		// error("a too large, ITMAX too small in routine gser");
		return NAN;
	}
}

// Returns the incomplete gamma function Q(a,x) evaluated by its continued fraction representation as gammcf.
// Also returns ln Γ(a) as gln.
double gcf( double *gammcf, double a, double x ){
	int i; 
	double an,b,c,d,del,h;
	
	double gln = gammln(a);
	b= x+1.0-a; 
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d; 
	for ( i = 0; i < ITMAX; i++ ) { // Iterate to convergence.
		an = -i*(i-a);
		b += 2.0; 
		d = an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c = b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d = 1.0/d;
		del = d*c; 
		h *= del;
		if (fabs(del-1.0) < DBL_EPSILON) break;
		//Set up for evaluating continued fraction by modified Lentz’s method (§5.2) with b0 = 0.
	} 
	if (i > ITMAX) error("a too large, ITMAX too small in gcf"); 
	*gammcf = exp(-x+a*log(x)-gln)*h;	//Put factors in front.
	return gln;
}

// Returns the value ln[Γ(xx)] for xx > 0.
double gammln( const double xx ) {
	// Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure accuracy is good enough.
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5}; 
	int j;
	y = x = xx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015; 
	for ( j = 0; j <= 5; j++ ) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

// Returns the value Γ(xx) for xx > 0
double gamm( const double xx ){
	return exp(gammln(xx));
}

// Returns x such that P.a;x/ D p for an argument p between 0 and 1.
double invgammp( const double p, const double a) {
	int j;
	double x,err,t,u,pp,a1=a-1; 
	double afac = 0;
	double lna1 = 0;
	double eps = 1.e-8; //Accuracy is the square of EPS.
	double gln = gammln(a);
	if (a <= 0.) return NAN; //error("a must be pos in invgammap"); 
	if (p >= 1.) return dmax(100.,a + 100.*sqrt(a)); 
	if (p <= 0.) return 0.0; 
	if (a > 1.) {
		lna1=log(a1);
		afac = exp(a1*(lna1-1.)-gln);
		pp = (p < 0.5)? p : 1. - p;
		t = sqrt(-2.*log(pp)); 
		x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
		if (p < 0.5) x = -x;
		x = dmax(1.e-3, a*pow(1.-1./(9.*a)-x/(3.*sqrt(a)),3));
	} 
	else { t = 1.0 - a*(0.253+a*0.12); 
		if (p < t) x = pow(p/t,1./a);
		else x = 1.-log(1.-(p-t)/(1.-t));
	} 
	for ( j = 0; j< 12; j++ ) {
		if (x <= 0.0) return 0.0; // x too small to compute accurately
		err = gammp(a,x) - p;
		if (a > 1.) t = afac*exp(-(x-a1)+a1*(log(x)-lna1)); 
		else t = exp(-x+a1*log(x)-gln);
		u = err/t;
		x -= (t = u/(1.-0.5*dmin(1.,u*((a-1.)/x - 1))));//	Halley’s method. 
		if (x <= 0.) x = 0.5*(x + t);	//Halve old value if x tries to go negative.
		if (fabs(t) < eps*x ) break;
	}
	return x;
}
