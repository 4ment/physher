//
//  beta.c
//  physher
//
//  Created by Mathieu Fourment on 15/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "beta.h"

#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "gamma.h"
#include "utils.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#define FPMIN 1.0e-30
#define EPS 3.0e-7

// Code from: Press, William H. et al. Numerical recipes 3rd edition: The art of scientific computing. 2007. Cambridge university press

//double betaiapprox(double a, double b, double x) {
////	Incomplete beta by quadrature. Returns Ix .a; b/. User should not call directly.
//	double xu,t,sum,ans;
//	double a1 = a-1.0, b1 = b-1.0, mu = a/(a+b);
//	double lnmu=log(mu),lnmuc=log(1.-mu);
//	t = sqrt(a*b/(SQR(a+b)*(a+b+1.0)));
//	if (x > a/(a+b)) { //Set how far to integrate into the tail:
//		if (x >= 1.0) return 1.0;
//		xu = dmin(1., dmax(mu + 10.*t, x + 5.0*t));
//	} else {
//		if (x <= 0.0) return 0.0;
//		xu = dmax(0., dmin(mu - 10.*t, x - 5.0*t));
//	}
//	sum = 0;
//	for (int j = 0; j < 18; j++) { // Gauss-Legendre.
//		t = x + (xu-x)*y[j];
//		sum += w[j]*exp(a1*(log(t)-lnmu)+b1*(log(1-t)-lnmuc));
//	}
//	ans = sum*(xu-x)*exp(a1*lnmu-gammln(a)+b1*lnmuc-gammln(b)+gammln(a+b));
//	return ans>0.0? 1.0-ans : -ans;
//}


double betacf(const double a, const double b, const double x) {
//	Evaluates continued fraction for incomplete beta function by modified Lentz’s method (􏱦5.2). User should not call directly.
//		272 Chapter 6. Special Functions
	int m2;
	double aa,c,d,del,h,qab,qam,qap;
	qab = a+b;
	qap = a+1.0;
	qam = a-1.0;
	c = 1.0;
	d = 1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d = 1.0/d;
	h = d;
	for (int m = 1; m < 10000; m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d = 1.0 + aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c = 1.0 + aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d = 1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d = 1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c = 1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d = 1.0/d;
		del = d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS) break;
	}
	return h;
}


double betai(const double a, const double b, const double x) {
	//Returns incomplete beta function Ix.a;b/ for positive a and b, and x between 0 and 1.
	if (a <= 0.0 || b <= 0.0){
		fprintf(stderr, "Bad a or b in routine betai");
		exit(1);
	}
	if (x < 0.0 || x > 1.0){
		fprintf(stderr, "Bad x in routine betai");
		exit(1);
	}
	if (x == 0.0 || x == 1.0) return x;
//	if (a > 3000 && b > 3000) return betaiapprox(a,b,x);
	double bt = exp(gammln(a+b) - gammln(a) - gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
	else return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double invbetai(double p, double a, double b) {
//	Inverse of incomplete beta function. Returns x such that Ix.a;b/ D p for argument p between 0 and 1.
//		const double EPS = 1.e-8;
	double pp,t,u,err,x,al,h,w,afac,a1=a-1.,b1=b-1.;
	int j;
	if (p <= 0.) return 0.;
	else if (p >= 1.) return 1.;
	else if (a >= 1. && b >= 1.) {
		pp = (p < 0.5)? p : 1. - p;
		t = sqrt(-2.*log(pp));
		x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
		if (p < 0.5) x = -x;
		al = (SQR(x)-3.)/6.;
		h = 2./(1./(2.*a-1.)+1./(2.*b-1.));
		w = (x*sqrt(al+h)/h)-(1./(2.*b-1)-1./(2.*a-1.))*(al+5./6.-2./(3.*h)); x = a/(a+b*exp(2.*w));
	} else {
		double lna = log(a/(a+b)), lnb = log(b/(a+b));
		t = exp(a*lna)/a;
		u = exp(b*lnb)/b;
		w = t + u;
		if (p < t/w) x = pow(a*w*p,1./a);
		else x = 1. - pow(b*w*(1.-p),1./b);
	}

	afac = -gammln(a) - gammln(b) + gammln(a+b);
	for (j=0;j<10;j++) {
		if (x == 0. || x == 1.) return x;
		err = betai(a,b,x) - p;
		t = exp(a1*log(x)+b1*log(1.-x) + afac);
		u = err/t; Halley:
		x -= (t = u/(1.-0.5*dmin(1.,u*(a1/x - b1/(1.-x)))));

		if (x <= 0.) x = 0.5*(x + t);
		if (x >= 1.) x = 0.5*(x + t + 1.); if (fabs(t) < EPS*x && j > 0) break;
	}
	return x;
}

double dbetaprime(double x, double alpha, double beta){
	 return pow(x, alpha - 1.0)*pow(1.0 + x, -alpha - beta)/gsl_sf_beta(alpha, beta);
}

double dlogbetaprime(double x, double alpha, double beta){
	return (alpha - 1.0)*log(x) - (alpha + beta)*log(1.0 + x) - gsl_sf_lnbeta(alpha, beta);
}
