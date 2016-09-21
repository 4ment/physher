/*
 *  neutralitytest.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 7/2/12.
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

#include "neutralitytest.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sequence.h"
#include "matrix.h"

static double calculate_pi( const Sequences *sequences ){
	double pi = 0.0;
	char *site1 = NULL;
	char *site2 = NULL;
	
	for ( int i = 0; i < sequences->size; i++ ) {
		for ( int j = i+1; j < sequences->size; j++ ) {
			site1 = sequences->seqs[i]->seq;
			site2 = sequences->seqs[j]->seq;
			for ( int k = 0; k < sequences->length; k++ ) {
				if ( *site1 != *site2 ) {
					pi++;
				}
				site1++;
				site2++;
			}
		}
	}
	double total = (sequences->size * (sequences->size-1) )*0.5; 
	return pi / total;
}

static int calculate_segragating_sites( const Sequences *sequences ){
	int i = 0;
	int j = 0;
	int S = 0;
	for ( int k = 0; k < sequences->length; k++ ) {
		for ( i = 0; i < sequences->size; i++ ) {
			for ( j = i+1; j < sequences->size; j++ ) {
				if ( sequences->seqs[i]->seq[k] != sequences->seqs[j]->seq[k] ) break;
			}
			if ( j != sequences->size ) {
				S++;
				break;
			}
		}
	}
	return S;
}

static int count_singleton( const Sequences *sequences ){
	int *counts = ivector(4);
	int nSingleton = 0;
	
	for (int i = 0; i < sequences->length; i++ ) {
		memset(counts, 0.0, 4*sizeof(int));
		for ( int j = 0; j < sequences->size; j++ ) {
			switch ( sequences->seqs[i]->seq[j] ) {
				case 'A':
					counts[0]++;
					break;
				case 'C':
					counts[1]++;
					break;
				case 'G':
					counts[2]++;
					break;
				case 'T':
					counts[3]++;
					break;
					
				default:
                    break;
			}
		}
		qsort( counts, 4, sizeof(int), qsort_desc_ivector );
		if ( counts [1] == 1 ) {
			nSingleton++;
		}
	}
	
	
	free(counts);
	return nSingleton;
}

// Null hypothesis E[pi] = theta and E[S] =a1*theta therefore pi=S/a1
// D is not normally distributed (almost)
double Tajima_D( const Sequences *sequences ){
	double a1 = 0.0;
	double a2 = 0.0;
	
	for ( int i = 1; i < sequences->size; i++ ) {
		a1 += 1.0/i;
		a2 += 1.0/(i*i);
	}
	double b1 = (sequences->size+1.0)/(3.0*(sequences->size-1));
	double b2 = 2.0*(sequences->size*sequences->size + sequences->size + 3)/(9.0*sequences->size*(sequences->size-1));
	double c1 = b1 - 1/a1;
	double c2 = b2 - (sequences->size+2)/(a1*sequences->size) + a2/(a1*a1);
	double e1 = c1/a1;
	double e2 = c2/(a1*a1+a2);
	double pi = calculate_pi(sequences);
	double S = calculate_segragating_sites(sequences);
	double D = ( pi - S/a1 )/(sqrt(e1*S + e2*S*(S-1)));
	
	fprintf(stderr, "Tajima's D\n\n");
	fprintf(stderr, "%d sequences, %d sites\n", sequences->size, sequences->length);
	fprintf(stderr, "pi:                %f\n", pi);
	fprintf(stderr, "Segragating sites: %d/%d \n", (int)S, sequences->length);
	fprintf(stderr, "Theta_hat[estimated from S]: %f\n", S/a1);
	fprintf(stderr, "Tajima's D: %f\n\n", D);
	
	fprintf(stderr, "a1: %f a2: %f b1: %f b2: %f\n", a1,a2,b1,b2);
	fprintf(stderr, "c1: %f c2: %f e1: %f e2: %f\n", c1,c2,e1,e2);
	
	return D;
}


// Theta = 4Nemu
// Ne: effective population size
// mu: mutation rate (per generation)
// Assumptions: n << Ne and infinitely many sites capable of varying
double Watterson_theta_estimator( const Sequences *sequences ){
	double K = calculate_segragating_sites(sequences);
	
	double a1 = 0.0;
	for ( int i = 1; i < sequences->size; i++ ) {
		a1 += 1.0/i;
	}
	double theta_hat = K/a1;
	return theta_hat;
}

// D*
double FuLi_Dstar( const Sequences *sequences ){
	double n = sequences->size;
	double nSingletons = count_singleton(sequences);
	double S = calculate_segragating_sites(sequences);
	
    double an = 0.0;
	double bn = 0.0;
    for (int i = 1; i < n; ++i){
		an += 1.0/i;
		bn += 1.0/(i*i);
	}
	
    double an_plus1 = an + 1.0/n;
	
    double cn = 2.0 * (n * an - 2.0 * (n - 1.0)) / ((n - 1.0) * (n - 2.0));
	
	
    double dn = cn + (n - 2.0) / (pow (n - 1.0, 2.0));
    dn += (2.0 / (n - 1.0)) * (1.5 - ((2.0 * an_plus1 - 3.0) / (n - 2.0)) -  1.0 / n);
	
	double vD = pow (n / (n - 1.0), 2.0) * bn;
    vD += an*an * dn;
    vD -= 2.0 * (n * an * (an + 1.0)) / (pow( (n - 1.0), 2.0) );
    vD /= an*an + bn;
	
	double uD = (n / (n - 1.0)) * (an - (n / (n - 1.0))) - vD;
	
	double Dstar = (n / (n - 1.0)) * S - an * nSingletons;
    Dstar /= sqrt( uD * S + vD * S * S );
	
	return Dstar;
}

// F*
double FuLi_Fstar( const Sequences *sequences ){
	double n = sequences->size;
	double nSingletons = count_singleton(sequences);
	double S = calculate_segragating_sites(sequences);
	double pi = calculate_pi(sequences);
	
    double an = 0.0;
	double bn = 0.0;
    for (int i = 1; i < n; ++i){
		an += 1.0/i;
		bn += 1.0/(i*i);
	}
	
    double an_plus1 = an + 1.0/n;
	
    double cn = 2.0 * (n * an - 2.0 * (n - 1.0)) / ((n - 1.0) * (n - 2.0));
	
	
    double dn = cn + (n - 2.0) / (pow (n - 1.0, 2.0));
    dn += (2.0 / (n - 1.0)) * (1.5 - ((2.0 * an_plus1 - 3.0) / (n - 2.0)) -  1.0 / n);
	
	double vF = (dn+2*(n*n+n+3)/(9.0*n*(n-1)) - 2.0/(n-1)*(4.0*bn-6.0+8.0/n))/(an*an+bn);
	double uF = (n/(n-1.0)+(n+1)/3/(n-1)-4/n/(n-1)+2*(n+1)/(n-1)/(n-1)*(an_plus1-2*n/(n+1)))/an-vF;
		
	double Fstar = pi - (n-1.0)/n*nSingletons;
    Fstar /= sqrt( uF * S + vF * S * S );
	
	return Fstar;
}


// test file for tajima D
//>a0
//ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA
//>a1
//AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA
//>a2
//AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA
//>a3
//AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA
//>a4
//AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA
//>a5
//AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA
//>a6
//AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA
//>a7
//AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA
//>a8
//AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA
//>a9
//AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA




