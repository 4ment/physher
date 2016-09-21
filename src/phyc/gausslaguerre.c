/*
 *  gausslaguerre.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 7/9/12.
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

#include <stdio.h>

#include <math.h>
#include "gamma.h"

#define EPS 3.0e-14
#define MAXIT 10

/* Given alf, the parameter α of the Laguerre polynomials, this routine returns arrays x[1..n] and w[1..n]
 * containing the abscissas and weights of the n-point Gauss-Laguerre quadrature formula. The smallest abscissa 
 * is returned in x[0], the largest in x[n-1].
 */
void gaulag( double  *x, double *w, int n, double alf ) {
	
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z1;
	double z = 0;
	
	double num = -exp(gammln(alf+n)-gammln((double)n));
	
	// Loop over the desired roots
	for ( i = 1; i <= n; i++ ) {
		// Initial guess for the smallest root
		if (i == 1) {
			z = (1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} 
		//Initial guess for the second root
		else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		}
		// Initial guess for the other roots
		else {
			ai = i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.0+3.5*ai))*(z-x[i-3])/(1.0+0.3*alf);
		}
		for ( its = 1; its <= MAXIT; its++ ) {
			p1 = 1.0;
			p2 = 0.0;
			for ( j = 1; j <= n; j++ ) {
				p3 = p2;
				p2 = p1;
				p1 = ((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			// p1 is now the desired Laguerre polynomial. We next compute pp,
			// its derivative, by a standard relation involving also p2, the polynomial of one lower order
			pp = (n*p1-(n+alf)*p2)/z;
			z1 = z;
			z = z1-p1/pp; // Newton’s formula.
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT){
			fprintf(stderr,"Too many iterations in gaulag\n");
		}
		x[i-1] = z; // Store the root and the weight. 
		w[i-1] = num/(pp*n*p2);
	}
}
