/*
 *  eigen.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/3/10.
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

#include "eigen.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#include "hessenberg.h"
#include "matrix.h"
#include "solve.h"

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
a[k][l]=h+s*(g-h*tau);

//#define JACOBI


#ifdef LAPACK_ENABLED
#include <lapacke.h>
#endif

static void hqr4(int n, int low, int hgh, double **h, double **zz, double *wr, double *wi);


static double complex _cdiv(double xr, double xi, double yr, double yi) {
    double r,d;
    double cdivr,cdivi;
    
    if (fabs(yr) > fabs(yi)) {
        r = yi/yr;
        d = yr + r*yi;
        cdivr = (xr + r*xi)/d;
        cdivi = (xi - r*xr)/d;
    }
    else {
        r = yr/yi;
        d = yi + r*yr;
        cdivr = (r*xr + xi)/d;
        cdivi = (r*xi - xr)/d;
    }
    double complex c = cdivr + cdivi*I;
    return  c;
}


EigenDecomposition * new_EigenDecomposition( const size_t dimension ){
	EigenDecomposition *ed = (EigenDecomposition *)malloc( sizeof(EigenDecomposition) );
	assert(ed);
	ed->dim = dimension;
	ed->eval  = dvector(dimension);
	ed->evali = dvector(dimension);
	ed->evec = dmatrix(dimension, dimension);
	ed->Invevec = dmatrix(dimension, dimension);
#ifdef LAPACK_ENABLED
    ed->isuppz = ivector(dimension);
    ed->M = dvector(dimension*dimension);
#endif
	return ed;
}

void free_EigenDecomposition( EigenDecomposition *eigencmp ){
	free(eigencmp->eval);
    free(eigencmp->evali);
	eigencmp->eval = NULL;
	free_dmatrix(eigencmp->evec, eigencmp->dim);
	free_dmatrix(eigencmp->Invevec, eigencmp->dim);
#ifdef LAPACK_ENABLED
    free(eigencmp->isuppz);
    free(eigencmp->M);
#endif
	free(eigencmp);
	eigencmp = NULL;
}

EigenDecomposition * clone_EigenDecomposition( EigenDecomposition *eigen ){
	EigenDecomposition *clone = (EigenDecomposition *)malloc( sizeof(EigenDecomposition) );
	assert(clone);
	clone->dim = eigen->dim;
	clone->evec    = clone_dmatrix(eigen->evec, eigen->dim, eigen->dim);
	clone->Invevec = clone_dmatrix(eigen->Invevec, eigen->dim, eigen->dim);
	clone->eval    = clone_dvector(eigen->eval, eigen->dim);
	clone->evali   = clone_dvector(eigen->evali, eigen->dim);
#ifdef LAPACK_ENABLED
    clone->isuppz = clone_ivector(eigen->isuppz, eigen->dim);
    clone->M      = clone_dvector(eigen->M, eigen->dim*eigen->dim);
#endif
	return clone;
}

/*
 lapack_int LAPACKE_dsyevr( int matrix_order, char jobz, char range, char uplo,
 lapack_int n, double* a, lapack_int lda, double vl,
 double vu, lapack_int il, lapack_int iu,
 double abstol, lapack_int* m, double* w, double* z,
 lapack_int ldz, lapack_int* isuppz );
 */
void EigenDecomposition_decompose( double **a, EigenDecomposition *eigen ){
	
	eigen->failed = false;
	
#ifdef LAPACK_ENABLED
    /* Negative abstol means using the default value */
    double abstol = -1.0;
    /* Set il, iu to compute NSELECT smallest eigenvalues */
    int il = 1;
    int iu = eigen->dim;
    
    int lda = eigen->dim;
    int ldz = eigen->dim;
    int n = eigen->dim;

    int info,m;
    double vl,vu;

    // should transpose it
    for( int i = 0; i < eigen->dim; i++ ){
        memcpy(&eigen->M[i*eigen->dim], &a[i], eigen->dim*sizeof(double));
    }

    // Solve eigenproblem
    // The matrix is row major but since the matrix is symmetric it is also column major
    // row major requires extra calculations (transpose)
    if( false ){
        info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'I', 'U', n, eigen->M, lda, vl, vu, il, iu, abstol, &m, eigen->eval, eigen->evec, ldz, eigen->isuppz);
    }
    else {
        double *H = clone_dvector(eigen->M, n*n);
        double *tau = dvector(n);
        printf("hess done...\n");
        
        // Reduce matrix to Hessenberg form
        info = LAPACKE_dgehrd(LAPACK_ROW_MAJOR, n, 1, n, H, lda, tau);
        printf("hess done...\n");
        
        // Generates the real orthogonal matrix Q
        info = LAPACKE_dorghr(LAPACK_COL_MAJOR, n, 1, n, eigen->M, n,tau);
        free(tau);
        printf("LAPACKE_orghr done...\n");
        
        // Computes all eigenvalues and the Schur factorization of a matrix reduced to Hessenberg form
        info = LAPACKE_dhseqr(LAPACK_COL_MAJOR, 'S', 'V', n, 1, n, H, n, eigen->eval, eigen->evali, eigen->M, n);
        printf("LAPACKE_dhseqr done...\n");
        
        // Computes right and left eigenvectors of an upper (quasi-) triangular matrix
        LAPACKE_dtrevc( LAPACK_COL_MAJOR, 'B', char howmny, lapack_logical* select, lapack_int n, H, lapack_int ldt, double* vl, lapack_int ldvl, double* vr, lapack_int ldvr, lapack_int mm, lapack_int* m);
        free(H);
    }
    
    /* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }
    printf("LAPACK\n");
    for ( int i = 0; i < eigen->dim; i++ ) {
        printf("%d %e \n", i, eigen->eval[i]);
    }
    exit(1);
    
    for ( int i = 0; i < n; i++ ) {
		memcpy(eigen->Invevec[i], eigen->evec[i], n*sizeof(double));
	}
    
    // LU decomoposition of a general matrix
    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, eigen->Invevec, lda, eigen->isuppz);
    if( info > 0 ) {
        printf( "The algorithm failed to compute LU decomoposition.\n" );
        exit( 1 );
    }
    
    
    // generate inverse of a matrix given its LU decomposition
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, eigen->Invevec, lda, eigen->isuppz);
    if( info > 0 ) {
        printf( "The algorithm failed to compute inverse from LU decomoposition.\n" );
        exit( 1 );
    }
    
#else
    
    int method = 1;
    
    size_t i;
    
    memset(eigen->eval, 0, sizeof(double)*eigen->dim);
    memset(eigen->evali, 0, sizeof(double)*eigen->dim);
    
    if(method == 0){
        double * d   = dvector(eigen->dim);
        double * b   = dvector(eigen->dim);
        int * index  = ivector(eigen->dim);
        int *order   = ivector(eigen->dim);
        
        
        for ( i = 0 ; i < eigen->dim; i++) {
            b[i]     = 1.0;
            d[i]     = 0.0;
            index[i] = 0;
            order[i] = 0;
        }	
    
    
        elmhes(a, order, eigen->dim );
        
        eltran( a, eigen->evec, order, eigen->dim);
        
        hqr4(eigen->dim, 1, eigen->dim, a, eigen->evec, eigen->eval, eigen->evali);
        for ( i = 0; i < eigen->dim; i++ ) {
            memcpy(eigen->Invevec[i], eigen->evec[i], eigen->dim*sizeof(double));
        }
        free(d);
        free(b);
        free(index);
        free(order);
    }
    else if( method == 1 ){
        double **H = clone_dmatrix(a, eigen->dim, eigen->dim);
        double *ort = dvector(eigen->dim);
        double **V  = dmatrix(eigen->dim, eigen->dim);
        
        orthes(eigen->dim, H, V, ort);
        
        int r = hqr2(eigen->dim, H, eigen->eval, eigen->evali, V, 10000);
        if( r == 0 ){
			eigen->failed = true;
//            print_dmatrix(stderr, (const double **)a, eigen->dim, eigen->dim, ' ');
//            error("Too many iterations in hqr2. Called from EigenDecomposition_decompose");
        }
        
        for ( i = 0; i < eigen->dim; i++ ) {
            memcpy(eigen->evec[i], V[i], eigen->dim*sizeof(double));
            memcpy(eigen->Invevec[i], V[i], eigen->dim*sizeof(double));
            free(H[i]);
            free(V[i]);
        }
        
        free(ort);
        free(H);
        free(V);
    }
	if(!eigen->failed){
		inverse(eigen->Invevec, eigen->dim);
	}
	
#endif
}


EigenDecomposition * eigen2( double **a, size_t dim ){
	
	//double **m = clone_dmatrix( a, dim, dim );
	
	int *order = ivector(dim);
	double ** vr = dmatrix(dim,dim);    //eigen vectors
	//double ** vi = dmatrix(dim,dim);  
	double * wr  = dvector(dim);  //eigen values
	double * wi  = dvector(dim);
	
	double * d  = dvector(dim);
	double * b  = dvector(dim);
	int * index  = ivector(dim);
	
	for ( size_t i = 0; i < dim; i++ ) {
		b[i] = 1.0;
		d[i] = 0.0;
		index[i] = 0;
	}
	
	EigenDecomposition *ed = NULL;
	ed = (EigenDecomposition *)malloc( sizeof(EigenDecomposition) );
	assert(ed);
	ed->dim = dim;

	
	/*fprintf(stdout, "A");
	print_matrix(a,4);*/

#ifndef JACOBI
	/*fprintf(stdout, "A");
	print_matrix(a,4);*/
	
    elmhes(a, order, dim );
	
	/*fprintf(stdout, "Elmhes");
	print_matrix(a,4);
	fprintf(stdout, "Order");
	for( int j = 0; j < dim; j++ ) fprintf(stdout, "%d  ",order[j]);
	fprintf(stdout, "\n");*/
	
	eltran( a, vr, order, dim);
	
	/*fprintf(stdout, "Eltran");
	print_matrix(a,4);
	fprintf(stdout, "Order");
	for( int j = 0; j < dim; j++ ) fprintf(stdout, "%d  ",order[j]);
	
	fprintf(stdout, "\n");
	
	fprintf(stdout, "Evec");
	print_matrix(vr,4);*/
	
	hqr4(dim, 1, dim, a, vr, wr, wi);
	
	/*print_matrix(vr,4);*/
	
#else
	jacobi( a, dim, wr, vr, order );
	
	fprintf(stdout, "Eigen values");
	for( size_t j = 0; j < dim; j++ ) fprintf(stdout, "%f  ",wr[j]);
	fprintf(stdout, "\n");
	print_matrix(vr,4);
#endif
	
	

	//fprintf(stdout, "\n");
	
	ed->evec = clone_dmatrix( vr, dim, dim );
	ed->eval = wr;
	
	inverse(vr, dim);
	
	/*fprintf(stdout, "Ivec\n");
	print_matrix(vr,4);*/
	
	ed->Invevec = vr;
	
	/*double **lambda = dmatrix(4,4);
	double **P = dmatrix(4,4);
	for ( i = 0; i < dim; i++ ){
		for ( int j = 0; j < dim; j++ ){
			lambda[i][j] = ( i == j ? ed->eval[i] : 0.);
		}
	}
	
	for ( i=0; i < ed->n; i++){
		for ( int j=0; j < ed->n; j++){
			P[i][j] = 0;
			for ( int  k = 0; k < ed->n; k++ )
				P[i][j] += ed->evec[i][k] * lambda[k][j];
		}
	}
	
	for ( i=0; i < ed->n; i++){
		for ( int j=0; j < ed->n; j++){
			lambda[i][j] = 0;
			for ( int  k = 0; k < ed->n; k++ )
				lambda[i][j] += P[i][k] * ed->Invevec[k][j];
		}
	}
	
	
	fprintf(stdout, "Id\n");
	print_matrix(lambda,4);*/
	
	free(order);
	free(index);
	free(wi);
	free(b);
	free(d);
    return ed; 
}


/* Computes all eigenvalues and eigenvectors of a real symmetric matrix a[0..n-1][0..n-1]. 
 * On output, elements of a above the diagonal are destroyed. d[0..n-1] returns the eigenvalues of a.
 * v[0..n-1][0..n-1] is a matrix whose columns contain, on output, the normalized eigenvectors of a.
 * nrot returns the number of Jacobi rotations that were required. */

void jacobi( double **a, int n, double *d, double **v, int *nrot){
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	b = dvector(n); 
	z = dvector(n); 
	//Initialize to the identity matrix.
	for ( ip = 0; ip < n; ip++) {
		for (iq= 0; iq < n; iq++) v[ip][iq]=0.0;
		v[ip][ip] = 1.0;
	}
	
	//Initialize b and d to the diagonal of a.
	for ( ip = 0; ip < n; ip++ ) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0; //This vector will accumulate terms of the form tapq as in equation (11.1.14).
	} 
	*nrot = 0;
	for ( i = 1; i <= 50; i++ ) {
		sm = 0.0;
		//Sum off-diagonal elements.
		for ( ip = 0; ip < n-1; ip++ ) {
			for ( iq = ip+1; iq < n; iq++ ) sm += fabs(a[ip][iq]);
		}
		//The normal return, which relies on quadratic convergence to machine underflow.
		if(sm == 0.0){
			free(z); 
			free(b);
			return;
		}
		//...on the first three sweeps.
		if (i < 4)
			tresh = 0.2*sm/(n*n); 
		else
			tresh=0.0; // ...thereafter.
		for ( ip = 0; ip < n-1; ip++ ){
			for ( iq = ip+1; iq < n; iq++ ) { 
				g = 100.0 * fabs(a[ip][iq]);
				//After four sweeps, skip the rotation if the off-diagonal element is small.
				if (i > 4 &&  (fabs(d[ip])+g) == fabs(d[ip]) && (fabs(d[iq])+g) == fabs(d[iq])) a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) { 
					h = d[iq]-d[ip];
					if ( (fabs(h)+g) == fabs(h)) t=(a[ip][iq])/h; //t = 1/(2θ)
					else { 
						theta = 0.5 * h/(a[ip][iq]);	//Equation (11.1.10).
						t = 1.0/(fabs(theta) + sqrt(1.0+theta*theta)); 
						if (theta < 0.0) t = -t;
					} 
					c = 1.0/sqrt(1+t*t);
					s=t*c; 
					tau = s/(1.0+c);
					h = t*a[ip][iq]; 
					z[ip] -= h; 
					z[iq] += h;
					d[ip] -= h; 
					d[iq] += h;
					a[ip][iq]=0.0; 
					
					for ( j = 0; j <= ip-1; j++ ) {     // Case of rotations 1 ≤ j < p.
						ROTATE(a,j,ip,j,iq)
					} 
					for (j = ip + 1; j <= iq-1; j++ ) { //Case of rotations p < j < q.
						ROTATE(a,ip,j,j,iq)
					}
					for ( j = iq + 1 ; j < n; j++ ) {   // Case of rotations q <j ≤n.
						ROTATE(a,ip,j,iq,j)
					} 
					for ( j = 0; j < n; j++ ) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				} 
			}
		}
		//Update d with the sum of tapq, and reinitialize z.
		for ( ip = 0; ip < n; ip++) {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] =0.0;
		}
	}
	error("Too many iterations in routine jacobi");
}

void normalize( double **a, size_t n ){
	double temp;
	size_t i,j;
	for ( i = 0; i < n ; i++) {
		temp = 0;
		for ( j = 0; j < n ; j++) temp += pow(a[i][j], 2);
		for ( j = 0; j < n ; j++) a[i][j] /= sqrt(temp);
	}
}

double complex mcdiv2(double ar, double ai, double br, double bi){
	double s, ars, ais, brs, bis;
	
	s = fabs(br) + fabs(bi);
	ars = ar/s;
	ais = ai/s;
	brs = br/s;
	bis = bi/s;
	s = brs * brs + bis * bis;
	double cr = (ars * brs + ais * bis)/s;
	double ci = (ais * brs - ars * bis)/s;
    double complex c = cr + ci*I;
    return c;
}

void hqr4(int n, int low, int hgh, double **h, double **zz, double *wr, double *wi) {
	int i, j, k, l=0, m, en, na, itn, its;
	double p=0, q=0, r=0, s=0, t, w, x=0, y, ra, sa, vi, vr, z=0, norm, tst1, tst2;
	bool notLast;
	
	
	norm = 0.0;
	k = 1;
	/* store isolated roots and compute matrix norm */
	for ( i = 0; i < n; i++ ) {
		for ( j = k - 1; j < n; j++ )	{
			norm += fabs(h[i][j]);
		}
		k = i + 1;
		if ( i + 1 < low || i + 1 > hgh ) {
			wr[i] = h[i][i];
			wi[i] = 0.0;
		}
	}
	en = hgh;
	t = 0.0;
	itn = n * 30;
	while (en >= low) {	/* search for next eigenvalues */
		its = 0;
		na = en - 1;
		while (en >= 1) {
			/* look for single small sub-diagonal element */
			bool fullLoop = true;
			for ( l = en; l > low; l-- ) {
				s = fabs(h[l - 2][l - 2]) + fabs(h[l - 1][l - 1]);
				if ( s == 0.0 ) {
					s = norm;
				}
				tst1 = s;
				tst2 = tst1 + fabs(h[l - 1][l - 2]);
				if ( tst2 == tst1 ) {
					fullLoop = false;
					break;
				}
			}
			if ( fullLoop ) {
				l = low;
			}
			
			x = h[en - 1][en - 1];	/* form shift */
			if ( l == en || l == na ) {
				break;
			}
			if ( itn == 0 ) {
				/* eigenvalues have not converged */
				error("Eigenvalues not converged");
			}
			y = h[na - 1][na - 1];
			w = h[en - 1][na - 1] * h[na - 1][en - 1];
			/* form exceptional shift */
			if ( its == 10 || its == 20 ) {
				t += x;
				for ( i = low - 1; i < en; i++ ) {
					h[i][i] -= x;
				}
				s = fabs(h[en - 1][na - 1]) + fabs(h[na - 1][en - 3]);
				x = 0.75 * s;
				y = x;
				w = -0.4375 * s * s;
			}
			its++;
			itn--;
			/* look for two consecutive small sub-diagonal elements */
			for ( m = en - 2; m >= l; m-- ) {
				z = h[m - 1][m - 1];
				r = x - z;
				s = y - z;
				p = (r * s - w) / h[m][m - 1] + h[m - 1][m];
				q = h[m][m] - z - r - s;
				r = h[m + 1][m];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if ( m == l ) {
					break;
				}
				tst1 = fabs(p) * (fabs(h[m - 2][m - 2]) + fabs(z) + fabs(h[m][m]));
				tst2 = tst1 + fabs(h[m - 1][m - 2]) * (fabs(q) + fabs(r));
				if (tst2 == tst1) {
					break;
				}
			}
			for ( i = m + 2; i <= en; i++ ) {
				h[i - 1][i - 3] = 0.0;
				if ( i != m + 2 ) {
					h[i - 1][i - 4] = 0.0;
				}
			}
			for ( k = m; k <= na; k++ )	{
				if ( k == na ) {
					notLast = false;
				}
				else {
					notLast = true;
				}
				if ( k != m ) {
					p = h[k - 1][k - 2];
					q = h[k][k - 2];
					r = 0.0;
					if ( notLast ) {
						r = h[k + 1][k - 2];
					}
					x = fabs(p) + fabs(q) + fabs(r);
					if ( x != 0.0 ) {
						p /= x;
						q /= x;
						r /= x;
					}
				}
				if ( x != 0.0 ) {
					if ( p < 0.0 ) {	/* sign */
						s = - sqrt(p * p + q * q + r * r);
					}
					else {
						s = sqrt(p * p + q * q + r * r);
					}
					if ( k != m ) {
						h[k - 1][k - 2] = -s * x;
					}
					else if ( l != m ) {
						h[k - 1][k - 2] = -h[k - 1][k - 2];
					}
					p += s;
					x = p / s;
					y = q / s;
					z = r / s;
					q /= p;
					r /= p;
					if ( !notLast ) {
						for ( j = k - 1; j < n; j++ ) {	/* row modification */
							p = h[k - 1][j] + q * h[k][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for (i = 0; i < j; i++) {	/* column modification */
							p = x * h[i][k - 1] + y * h[i][k];
							h[i][k - 1] -= p;
							h[i][k] -= p * q;
						}
						/* accumulate transformations */
						for (i = low - 1; i < hgh; i++) {
							p = x * zz[i][k - 1] + y * zz[i][k];
							zz[i][k - 1] -= p;
							zz[i][k] -= p * q;
						}
					}
					else
					{
						for ( j = k - 1; j < n; j++ ) {	/* row modification */
							p = h[k - 1][j] + q * h[k][j] + r * h[k + 1][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
							h[k + 1][j] -= p * z;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for ( i = 0; i < j; i++) {	/* column modification */
							p = x * h[i][k - 1] + y * h[i][k] + z * h[i][k + 1];
							h[i][k - 1] -= p;
							h[i][k] -= p * q;
							h[i][k + 1] -= p * r;
						}
						/* accumulate transformations */
						for ( i = low - 1; i < hgh; i++ ) {
							p = x * zz[i][k - 1] + y * zz[i][k] +
							z * zz[i][k + 1];
							zz[i][k - 1] -= p;
							zz[i][k] -= p * q;
							zz[i][k + 1] -= p * r;
						}
					}
				}
			}				 /* for k */
		}					 /* while infinite loop */
		if ( l == en ) {				 /* one root found */
			h[en - 1][en - 1] = x + t;
			wr[en - 1] = h[en - 1][en - 1];
			wi[en - 1] = 0.0;
			en = na;
			continue;
		}
		y = h[na - 1][na - 1];
		w = h[en - 1][na - 1] * h[na - 1][en - 1];
		p = (y - x) / 2.0;
		q = p * p + w;
		z = sqrt(fabs(q));
		h[en - 1][en - 1] = x + t;
		x = h[en - 1][en - 1];
		h[na - 1][na - 1] = y + t;
		if (q >= 0.0) {	 /* real pair */
			if (p < 0.0) {	/* sign */
				z = p - fabs(z);
			}
			else {
				z = p + fabs(z);
			}
			wr[na - 1] = x + z;
			wr[en - 1] = wr[na - 1];
			if (z != 0.0) {
				wr[en - 1] = x - w / z;
			}
			wi[na - 1] = 0.0;
			wi[en - 1] = 0.0;
			x = h[en - 1][na - 1];
			s = fabs(x) + fabs(z);
			p = x / s;
			q = z / s;
			r = sqrt(p * p + q * q);
			p /= r;
			q /= r;
			for (j = na - 1; j < n; j++)
			{	/* row modification */
				z = h[na - 1][j];
				h[na - 1][j] = q * z + p * h[en - 1][j];
				h[en - 1][j] = q * h[en - 1][j] - p * z;
			}
			for (i = 0; i < en; i++) {	/* column modification */
				z = h[i][na - 1];
				h[i][na - 1] = q * z + p * h[i][en - 1];
				h[i][en - 1] = q * h[i][en - 1] - p * z;
			}
			/* accumulate transformations */
			for (i = low - 1; i < hgh; i++) {
				z = zz[i][na - 1];
				zz[i][na - 1] = q * z + p * zz[i][en - 1];
				zz[i][en - 1] = q * zz[i][en - 1] - p * z;
			}
		}
		else {	/* complex pair */
			wr[na - 1] = x + p;
			wr[en - 1] = x + p;
			wi[na - 1] = z;
			wi[en - 1] = -z;
		}
		en -= 2;
	} /* while en >= low */
	/* backsubstitute to find vectors of upper triangular form */
	if (norm != 0.0) {
		for (en = n; en >= 1; en--) {
			p = wr[en - 1];
			q = wi[en - 1];
			na = en - 1;
			if (q == 0.0) {/* real vector */
				m = en;
				h[en - 1][en - 1] = 1.0;
				if (na != 0) {
					for (i = en - 2; i >= 0; i--) {
						w = h[i][i] - p;
						r = 0.0;
						for (j = m - 1; j < en; j++) {
							r += h[i][j] * h[j][en - 1];
						}
						if (wi[i] < 0.0) {
							z = w;
							s = r;
						}
						else {
							m = i + 1;
							if (wi[i] == 0.0) {
								t = w;
								if (t == 0.0) {
									tst1 = norm;
									t = tst1;
									do {
										t = 0.01 * t;
										tst2 = norm + t;
									}
									while (tst2 > tst1);
								}
								h[i][en - 1] = -(r / t);
							}
							else {	/* solve real equations */
								x = h[i][i + 1];
								y = h[i + 1][i];
								q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
								t = (x * s - z * r) / q;
								h[i][en - 1] = t;
								if (fabs(x) > fabs(z))
									h[i + 1][en - 1] = (-r - w * t) / x;
								else
									h[i + 1][en - 1] = (-s - y * t) / z;
							}
							/* overflow control */
							t = fabs(h[i][en - 1]);
							if (t != 0.0) {
								tst1 = t;
								tst2 = tst1 + 1.0 / tst1;
								if (tst2 <= tst1) {
									for (j = i; j < en; j++) {
										h[j][en - 1] /= t;
									}
								}
							}
						}
					}
				}
			}
			else if (q > 0.0) {
				m = na;
				if ( fabs(h[en - 1][na - 1]) > fabs(h[na - 1][en - 1]) ) {
					h[na - 1][na - 1] = q / h[en - 1][na - 1];
					h[na - 1][en - 1] = (p - h[en - 1][en - 1])/h[en - 1][na - 1];
				}
				else {
					double complex z1 = mcdiv2(0.0, -h[na - 1][en - 1], h[na - 1][na - 1] - p, q);
					h[na - 1][na - 1] = creal(z1);
					h[na - 1][en - 1] = cimag(z1);
				}
				h[en - 1][na - 1] = 0.0;
				h[en - 1][en - 1] = 1.0;
				if (en != 2) {
					for (i = en - 3; i >= 0; i--) {
						w = h[i][i] - p;
						ra = 0.0;
						sa = 0.0;
						for (j = m - 1; j < en; j++) {
							ra += h[i][j] * h[j][na - 1];
							sa += h[i][j] * h[j][en - 1];
						}
						if (wi[i] < 0.0) {
							z = w;
							r = ra;
							s = sa;
						}
						else {
							m = i + 1;
							if (wi[i] == 0.0) {
								double complex z = mcdiv2(-ra, -sa, w, q);
								h[i][na - 1] = creal(z);
								h[i][en - 1] = cimag(z);
							}
							else {	/* solve complex equations */
								x = h[i][i + 1];
								y = h[i + 1][i];
								vr = (wr[i] - p) * (wr[i] - p);
								vr = vr + wi[i] * wi[i] - q * q;
								vi = (wr[i] - p) * 2.0 * q;
								if (vr == 0.0 && vi == 0.0) {
									tst1 = norm * (fabs(w) + fabs(q) + fabs(x) +
												   fabs(y) + fabs(z));
									vr = tst1;
									do {
										vr = 0.01 * vr;
										tst2 = tst1 + vr;
									}
									while (tst2 > tst1);
								}
								double complex z1 = mcdiv2(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
								h[i][na - 1] = creal(z1);
								h[i][en - 1] = cimag(z1);
								if (fabs(x) > fabs(z) + fabs(q)) {
									h[i + 1]
									[na - 1] = (q * h[i][en - 1] -
												w * h[i][na - 1] - ra) / x;
									h[i + 1][en - 1] = (-sa - w * h[i][en - 1] -
														q * h[i][na - 1]) / x;
								}
								else {
									double complex z1 = mcdiv2(-r - y * h[i][na - 1], -s - y * h[i][en - 1], z, q);
									h[i + 1][na - 1] = creal(z1);
									h[i + 1][en - 1] = cimag(z1);
								}
							}
							/* overflow control */
							t = (fabs(h[i][na - 1]) > fabs(h[i][en - 1])) ?
							fabs(h[i][na - 1]) : fabs(h[i][en - 1]);
							if (t != 0.0) {
								tst1 = t;
								tst2 = tst1 + 1.0 / tst1;
								if (tst2 <= tst1) {
									for (j = i; j < en; j++) {
										h[j][na - 1] /= t;
										h[j][en - 1] /= t;
									}
								}
							}
						}
					}
				}
			}
		}
		/* end back substitution. vectors of isolated roots */
		for ( i = 0; i < n; i++ ) {
			if (i + 1 < low || i + 1 > hgh ) {
				for ( j = i; j < n; j++ ) {
					zz[i][j] = h[i][j];
				}
			}
		}
		/* multiply by transformation matrix to give vectors of
		 * original full matrix. */
		for ( j = n - 1; j >= low - 1; j-- ) {
			m = ((j + 1) < hgh) ? (j + 1) : hgh; /* min */
			for ( i = low - 1; i < hgh; i++ ) {
				z = 0.0;
				for ( k = low - 1; k < m; k++ ) {
					z += zz[i][k] * h[k][j];
				}
				zz[i][j] = z;
			}
		}
	}
}


// Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n] by qr method.
//  This is derived from the Algol procedure hqr2,
//  by Martin and Wilkinson, Handbook for Auto. Comp.,
//  Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.
int hqr2( int N, double **H, double *d, double *e, double **V, int maxIterations ) {
    
    // Initialize
    
    int nn = N;
    int n = nn-1;
    int low = 0;
    int high = nn-1;
    double eps = pow(2.0,-52.0);
    double exshift = 0.0;
    double p=0,q=0,r=0,s=0,z=0,t,w,x,y;
    
    // Store roots isolated by balanc and compute matrix norm
    
    double norm = 0.0;
    for (int i = 0; i < nn; i++) {
        if (i < low | i > high) {
            d[i] = H[i][i];
            e[i] = 0.0;
        }
        for (int j = imax(i-1,0); j < nn; j++) {
            norm = norm + fabs(H[i][j]);
        }
    }
    
    // Outer loop over eigenvalue index
    
    int iter = 0;
    while (n >= low) {
        
        // Look for single small sub-diagonal element
        
        int l = n;
        while (l > low) {
            s = fabs(H[l-1][l-1]) + fabs(H[l][l]);
            if (s == 0.0) {
                s = norm;
            }
            if (fabs(H[l][l-1]) < eps * s) {
                break;
            }
            l--;
        }
        
        // Check for convergence
        // One root found
        
        if (l == n) {
            H[n][n] = H[n][n] + exshift;
            d[n] = H[n][n];
            e[n] = 0.0;
            n--;
            iter = 0;
            
            // Two roots found
            
        } else if (l == n-1) {
            w = H[n][n-1] * H[n-1][n];
            p = (H[n-1][n-1] - H[n][n]) / 2.0;
            q = p * p + w;
            z = sqrt(fabs(q));
            H[n][n] = H[n][n] + exshift;
            H[n-1][n-1] = H[n-1][n-1] + exshift;
            x = H[n][n];
            
            // Real pair
            
            if (q >= 0) {
                if (p >= 0) {
                    z = p + z;
                } else {
                    z = p - z;
                }
                d[n-1] = x + z;
                d[n] = d[n-1];
                if (z != 0.0) {
                    d[n] = x - w / z;
                }
                e[n-1] = 0.0;
                e[n] = 0.0;
                x = H[n][n-1];
                s = fabs(x) + fabs(z);
                p = x / s;
                q = z / s;
                r = sqrt(p * p+q * q);
                p = p / r;
                q = q / r;
                
                // Row modification
                
                for (int j = n-1; j < nn; j++) {
                    z = H[n-1][j];
                    H[n-1][j] = q * z + p * H[n][j];
                    H[n][j] = q * H[n][j] - p * z;
                }
                
                // Column modification
                
                for (int i = 0; i <= n; i++) {
                    z = H[i][n-1];
                    H[i][n-1] = q * z + p * H[i][n];
                    H[i][n] = q * H[i][n] - p * z;
                }
                
                // Accumulate transformations
                
                for (int i = low; i <= high; i++) {
                    z = V[i][n-1];
                    V[i][n-1] = q * z + p * V[i][n];
                    V[i][n] = q * V[i][n] - p * z;
                }
                
                // Complex pair
                
            } else {
                d[n-1] = x + p;
                d[n] = x + p;
                e[n-1] = z;
                e[n] = -z;
            }
            n = n - 2;
            iter = 0;
            
            // No convergence yet
            
        } else {
            
            // Form shift
            
            x = H[n][n];
            y = 0.0;
            w = 0.0;
            if (l < n) {
                y = H[n-1][n-1];
                w = H[n][n-1] * H[n-1][n];
            }
            
            // Wilkinson's original ad hoc shift
            
            if (iter == 10) {
                exshift += x;
                for (int i = low; i <= n; i++) {
                    H[i][i] -= x;
                }
                s = fabs(H[n][n-1]) + fabs(H[n-1][n-2]);
                x = y = 0.75 * s;
                w = -0.4375 * s * s;
            }
            
            // MATLAB's new ad hoc shift
            
            if (iter == 30) {
                s = (y - x) / 2.0;
                s = s * s + w;
                if (s > 0) {
                    s = sqrt(s);
                    if (y < x) {
                        s = -s;
                    }
                    s = x - w / ((y - x) / 2.0 + s);
                    for (int i = low; i <= n; i++) {
                        H[i][i] -= s;
                    }
                    exshift += s;
                    x = y = w = 0.964;
                }
            }
            
            iter = iter + 1;   // (Could check iteration count here.)
            if(iter > maxIterations){
                return 0;
                //error("Too many iterations in hqr2");
            }
            
            // Look for two consecutive small sub-diagonal elements
            
            int m = n-2;
            while (m >= l) {
                z = H[m][m];
                r = x - z;
                s = y - z;
                p = (r * s - w) / H[m+1][m] + H[m][m+1];
                q = H[m+1][m+1] - z - r - s;
                r = H[m+2][m+1];
                s = fabs(p) + fabs(q) + fabs(r);
                p = p / s;
                q = q / s;
                r = r / s;
                if (m == l) {
                    break;
                }
                if (fabs(H[m][m-1]) * (fabs(q) + fabs(r)) <
                    eps * (fabs(p) * (fabs(H[m-1][m-1]) + fabs(z) +
                                          fabs(H[m+1][m+1])))) {
                    break;
                }
                m--;
            }
            
            for (int i = m+2; i <= n; i++) {
                H[i][i-2] = 0.0;
                if (i > m+2) {
                    H[i][i-3] = 0.0;
                }
            }
            
            // Double QR step involving rows l:n and columns m:n
            
            for (int k = m; k <= n-1; k++) {
                bool notlast = (k != n-1);
                if (k != m) {
                    p = H[k][k-1];
                    q = H[k+1][k-1];
                    r = (notlast ? H[k+2][k-1] : 0.0);
                    x = fabs(p) + fabs(q) + fabs(r);
                    if (x != 0.0) {
                        p = p / x;
                        q = q / x;
                        r = r / x;
                    }
                }
                if (x == 0.0) {
                    break;
                }
                s = sqrt(p * p + q * q + r * r);
                if (p < 0) {
                    s = -s;
                }
                if (s != 0) {
                    if (k != m) {
                        H[k][k-1] = -s * x;
                    } else if (l != m) {
                        H[k][k-1] = -H[k][k-1];
                    }
                    p = p + s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q = q / p;
                    r = r / p;
                    
                    // Row modification
                    
                    for (int j = k; j < nn; j++) {
                        p = H[k][j] + q * H[k+1][j];
                        if (notlast) {
                            p = p + r * H[k+2][j];
                            H[k+2][j] = H[k+2][j] - p * z;
                        }
                        H[k][j] = H[k][j] - p * x;
                        H[k+1][j] = H[k+1][j] - p * y;
                    }
                    
                    // Column modification
                    
                    for (int i = 0; i <= imin(n,k+3); i++) {
                        p = x * H[i][k] + y * H[i][k+1];
                        if (notlast) {
                            p = p + z * H[i][k+2];
                            H[i][k+2] = H[i][k+2] - p * r;
                        }
                        H[i][k] = H[i][k] - p;
                        H[i][k+1] = H[i][k+1] - p * q;
                    }
                    
                    // Accumulate transformations
                    
                    for (int i = low; i <= high; i++) {
                        p = x * V[i][k] + y * V[i][k+1];
                        if (notlast) {
                            p = p + z * V[i][k+2];
                            V[i][k+2] = V[i][k+2] - p * r;
                        }
                        V[i][k] = V[i][k] - p;
                        V[i][k+1] = V[i][k+1] - p * q;
                    }
                }  // (s != 0)
            }  // k loop
        }  // check convergence
    }  // while (n >= low)
    
    // Backsubstitute to find vectors of upper triangular form
    
    if (norm == 0.0) {
        return -1;
    }
    
    for (n = nn-1; n >= 0; n--) {
        p = d[n];
        q = e[n];
        
        // Real vector
        
        if (q == 0) {
            int l = n;
            H[n][n] = 1.0;
            for (int i = n-1; i >= 0; i--) {
                w = H[i][i] - p;
                r = 0.0;
                for (int j = l; j <= n; j++) {
                    r = r + H[i][j] * H[j][n];
                }
                if (e[i] < 0.0) {
                    z = w;
                    s = r;
                } else {
                    l = i;
                    if (e[i] == 0.0) {
                        if (w != 0.0) {
                            H[i][n] = -r / w;
                        } else {
                            H[i][n] = -r / (eps * norm);
                        }
                        
                        // Solve real equations
                        
                    } else {
                        x = H[i][i+1];
                        y = H[i+1][i];
                        q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                        t = (x * s - z * r) / q;
                        H[i][n] = t;
                        if (fabs(x) > fabs(z)) {
                            H[i+1][n] = (-r - w * t) / x;
                        } else {
                            H[i+1][n] = (-s - y * t) / z;
                        }
                    }
                    
                    // Overflow control
                    
                    t = fabs(H[i][n]);
                    if ((eps * t) * t > 1) {
                        for (int j = i; j <= n; j++) {
                            H[j][n] = H[j][n] / t;
                        }
                    }
                }
            }
            
            // Complex vector
            
        } else if (q < 0) {
            int l = n-1;
            // Last vector component imaginary so matrix is triangular
            
            if (fabs(H[n][n-1]) > fabs(H[n-1][n])) {
                H[n-1][n-1] = q / H[n][n-1];
                H[n-1][n] = -(H[n][n] - p) / H[n][n-1];
            } else {
                double complex c = _cdiv(0.0,-H[n-1][n],H[n-1][n-1]-p,q);
                H[n-1][n-1] = creal(c);
                H[n-1][n] = cimag(c);
            }
            H[n][n-1] = 0.0;
            H[n][n] = 1.0;
            for (int i = n-2; i >= 0; i--) {
                double ra,sa,vr,vi;
                ra = 0.0;
                sa = 0.0;
                for (int j = l; j <= n; j++) {
                    ra = ra + H[i][j] * H[j][n-1];
                    sa = sa + H[i][j] * H[j][n];
                }
                w = H[i][i] - p;
                
                if (e[i] < 0.0) {
                    z = w;
                    r = ra;
                    s = sa;
                } else {
                    l = i;
                    if (e[i] == 0) {
                        double complex c = _cdiv(-ra,-sa,w,q);
                        H[i][n-1] = creal(c);
                        H[i][n] = cimag(c);
                    } else {
                        
                        // Solve complex equations
                        
                        x = H[i][i+1];
                        y = H[i+1][i];
                        vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                        vi = (d[i] - p) * 2.0 * q;
                        if (vr == 0.0 & vi == 0.0) {
                            vr = eps * norm * (fabs(w) + fabs(q) +
                                               fabs(x) + fabs(y) + fabs(z));
                        }
                        double complex c = _cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
                        H[i][n-1] = creal(c);
                        H[i][n] = cimag(c);
                        if (fabs(x) > (fabs(z) + fabs(q))) {
                            H[i+1][n-1] = (-ra - w * H[i][n-1] + q * H[i][n]) / x;
                            H[i+1][n] = (-sa - w * H[i][n] - q * H[i][n-1]) / x;
                        } else {
                            double complex c = _cdiv(-r-y*H[i][n-1],-s-y*H[i][n],z,q);
                            H[i+1][n-1] = creal(c);
                            H[i+1][n] = cimag(c);
                        }
                    }
                    
                    // Overflow control
                    
                    t = imax(fabs(H[i][n-1]),fabs(H[i][n]));
                    if ((eps * t) * t > 1) {
                        for (int j = i; j <= n; j++) {
                            H[j][n-1] = H[j][n-1] / t;
                            H[j][n] = H[j][n] / t;
                        }
                    }
                }
            }
        }
    }
    
    // Vectors of isolated roots
    
    for (int i = 0; i < nn; i++) {
        if (i < low | i > high) {
            for (int j = i; j < nn; j++) {
                V[i][j] = H[i][j];
            }
        }
    }
    
    // Back transformation to get eigenvectors of original matrix
    
    for (int j = nn-1; j >= low; j--) {
        for (int i = low; i <= high; i++) {
            z = 0.0;
            for (int k = low; k <= imin(j,high); k++) {
                z = z + V[i][k] * H[k][j];
            }
            V[i][j] = z;
        }
    }
    return 1;
}

//// Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n].
//// On input a can be exactly as output from elmhes §11.5; on output it is destroyed.
//// The real and imaginary parts of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
//void hqr2(double **a, int n, double *wr, double *wi, double **vr, double **vi, int low, int hi ){
//    int nn,m,l,k,j,its,i,mmin;
//    double z,y,x,w,vv,uu,t,s,r,q,p,anorm;
//    s= r = p = q = z = anorm = 0.0;
//    
//    double ra;
//    complex v;
//    double sa,eps = 0.0;
//    
//    k = 1;
//    // Compute matrix norm for possible use in locating single small subdiagonal element.
//    for ( i = 0; i < n; i++ ){
//        for ( j = imax(i-1,0); j < n; j++ ){
//            anorm += fabs(a[i][j]);
//        }
//        k = i + 1;
//        if (i < low || i + 1 > hi) {
//            wr[i] = a[i][i];
//            wi[i] = 0.0;
//        }
//    }
//    //nn = n-1;
//    nn = hi-1;
//    t=0.0;
//    // Gets changed only by an exceptional shift.
//    // Begin search for next eigenvalue.
//    while(nn >= low) {
//        its=0;
//        do{
//            // Begin iteration: look for single small subdiagonal element.
//            for( l = nn; l > low; l-- ) {
//                s = fabs(a[l-1][l-1]) + fabs(a[l][l]);
//                if (s == 0.0) s = anorm;
//                if ( fabs(a[l][l-1] + s) == s) break;
//            }
//            
//            x = a[nn][nn];
//            
//            //One root found.
//            if (l == nn) {
//                wr[nn] = x+t;
//                a[nn][nn] = x + t; // for eigenvector
//                wi[nn--] = 0.0;
//            }
//            else {
//                y = a[nn-1][nn-1];
//                w = a[nn][nn-1]*a[nn-1][nn];
//                
//                //Two roots found...
//                if (l == nn-1 ) {
//                    p = 0.5*(y-x);
//                    q = p*p+w;
//                    z = sqrt(fabs(q));
//                    x += t;
//                    
//                    a[nn][nn] = x; // for eigenvector
//                    a[nn-1][nn-1] = y + t;
//                    
//                    //	...a real pair.
//                    if (q >= 0.0) {
//                        z =  p + SIGN(z,p);
//                        wr[nn-1] = wr[nn] = x+z;
//                        if (z) wr[nn] = x-w/z;
//                        wi[nn-1] = wi[nn]=0.0;
//                        
//                        // for egeinvector
//                        x = a[nn][nn-1];
//                        s = fabs(x) + fabs(z);
//                        p = x / s;
//                        q = z / s;
//                        r = sqrt(p*p+q*q);
//                        p /= r;
//                        q /= r;
//                        for (j = nn - 1; j < n; j++) {
//                            z = a[nn-1][j];
//                            a[nn-1][j] = q * z + p *
//                            a[nn][j];
//                            a[nn][j] = q * a[nn][j] - p*z;
//                        }
//                        for (i = 0; i <= nn; i++) {
//                            z = a[i][nn-1];
//                            a[i][nn-1] = q * z + p * a[i][nn];
//                            a[i][nn] = q * a[i][nn] - p*z;
//                        }
//                        for (i = low; i < hi; i++) {
//                            z = vr[i][nn-1];
//                            vr[i][nn-1] = q*z + p*vr[i][nn];
//                            vr[i][nn] = q*vr[i][nn] - p*z;
//                        }
//                        
//                    }
//                    //...a complex pair.
//                    else{
//                        wr[nn-1] = wr[nn]=x+p;
//                        wi[nn-1] = -(wi[nn]=z);
//                    }
//                    nn -= 2;
//                }
//                // No roots found. Continue iteration.
//                else{
//                    if (its == 30) error("Too many iterations in hqr");
//                    // Form exceptional shift.
//                    if (its == 10 || its == 20) {
//                        t += x;
//                        for ( i=low; i<=nn; i++ ) a[i][i] -= x;
//                        s = fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
//                        y = x = 0.75*s;
//                        w = -0.4375*s*s;
//                    }
//                    ++its;
//                    //Form shift and then look for 2 consecutive small sub-diagonal elements.
//                    for( m=(nn-2); m>=l; m-- ) {
//                        z = a[m][m];
//                        r = x-z;
//                        s = y-z;
//                        p = (r*s-w)/a[m+1][m]+a[m][m+1]; //Equation (11.6.23).
//                        q = a[m+1][m+1]-z-r-s;
//                        r = a[m+2][m+1];
//                        
//                        
//                        // Scale to prevent overflow or underflow.
//                        s=fabs(p)+fabs(q)+fabs(r);
//                        
//                        p /= s;
//                        q /= s;
//                        r /= s;
//                        if (m == l) break;
//                        uu = fabs(a[m][m-1])*(fabs(q)+fabs(r));
//                        vv = fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1])); //Equation (11.6.26).
//                        if( (double)(uu+vv) == vv) break;
//                    }
//                    for ( i = m+2; i <= nn; i++ ) {
//                        a[i][i-2]=0.0;
//                        if (i != (m+2)) a[i][i-3]=0.0;
//                    }
//                    for ( k = m; k < nn; k++ ) {
//                        // Double QR step on rows l to nn and columns m to nn.
//                        if (k != m) {
//                            p = a[k][k-1];	//Begin setup of Householder vector
//                            q = a[k+1][k-1];
//                            r=0.0;
//                            if (k != (nn-1)) r=a[k+2][k-1];
//                            if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
//                                //Scale to prevent overflow or underflow.
//                                p /= x;
//                                q /= x;
//                                r /= x;
//                            }
//                        }
//                        if ((s = SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
//                            if (k == m) {
//                                if (l != m)
//                                    a[k][k-1] = -a[k][k-1];
//                            }
//                            else a[k][k-1] = -s*x;
//                            p += s;     //Equations (11.6.24).
//                            x = p/s;
//                            y = q/s;
//                            z = r/s;
//                            q /= p;
//                            r /= p;
//                            
//                            // Row modification.
//                            for ( j = k; j < n-1; j++ ) {
//                                p = a[k][j]+q*a[k+1][j];
//                                if (k != (nn-1)) {
//                                    p += r*a[k+2][j];
//                                    a[k+2][j] -= p*z;
//                                }
//                                a[k+1][j] -= p*y;
//                                a[k][j]   -= p*x;
//                            }
//                            mmin = nn < k+3 ? nn : k+3;
//                            
//                            // Column modification.
//                            //for ( i = l; i <= mmin; i++ ) { // for egeinvector
//                            for ( i = 0; i <= mmin; i++ ) {
//                                p = x*a[i][k]+y*a[i][k+1];
//                                if (k != (nn-1)) {
//                                    p += z*a[i][k+2];
//                                    a[i][k+2] -= p*r;
//                                }
//                                a[i][k+1] -= p*q;
//                                a[i][k]   -= p;
//                            }
//                            
//                            // for egeinvector
//                            for (i = low; i < hi; i++) {
//                                p = x * vr[i][k] + y * vr[i][k+1];
//                                if (k != nn - 1) {
//                                    p += z * vr[i][k+2];
//                                    vr[i][k+2] -= p*r;
//                                }
//                                vr[k+1][n] -= p*q;
//                                vr[i][k] -= p;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        while( 1 < nn-1);
//    }
//    
//    if (anorm != 0) {
//        /* back substitute to find vectors of upper triangular form */
//        for ( nn = n-1; nn >= 0; nn--) {
//            p = wr[nn];
//            if ((q = wi[nn]) < 0) {            /* complex vector */
//                m = nn - 1;
//                if ( fabs(a[nn][nn-1]) > fabs(a[nn-1][nn]) ) {
//                    a[nn-1][nn-1] = q / a[nn][nn-1];
//                    a[nn-1][nn] = (p - a[nn][nn]) / a[nn][nn-1];
//                }
//                else {
//                    v = cdiv(compl(0.0,-a[nn-1][nn]),
//                             compl(a[nn-1][nn-1]-p,q));
//                    a[nn-1][nn-1] = v.re;
//                    a[nn-1][nn] = v.im;
//                }
//                a[nn][nn-1] = 0;
//                a[nn][nn] = 1;
//                for (i = nn - 2; i >= 0; i--) {
//                    w = a[i][i] - p;
//                    ra = 0;
//                    sa = a[i][nn];
//                    for (j = m; j < nn; j++) {
//                        ra += a[i][j] * a[j][nn-1];
//                        sa += a[i][j] * a[j][nn];
//                    }
//                    if (wi[i] < 0) {
//                        z = w;
//                        r = ra;
//                        s = sa;
//                    }
//                    else {
//                        m = i;
//                        if (wi[i] == 0) {
//                            v = cdiv(compl(-ra,-sa),compl(w,q));
//                            a[i][nn-1] = v.re;
//                            a[i][nn] = v.im;
//                        }
//                        else {                      /* solve complex equations */
//                            x = a[i][i+1];
//                            y = a[i+1][i];
//                            v.re = (wr[i]- p)*(wr[i]-p) + wi[i]*wi[i] - q*q;
//                            v.im = (wr[i] - p)*2*q;
//                            if ((fabs(v.re) + fabs(v.im)) == 0) {
//                                v.re = eps * anorm * (fabs(w) +
//                                                      fabs(q) + fabs(x) + fabs(y) + fabs(z));
//                            }
//                            v = cdiv(compl(x*r-z*ra+q*sa,x*s-z*sa-q*ra),v);
//                            a[i][nn-1] = v.re;
//                            a[i][nn] = v.im;
//                            if (fabs(x) > fabs(z) + fabs(q)) {
//                                a[i+1][nn-1] =  (-ra - w * a[i][nn-1] + q * a[i][nn]) / x;
//                                a[i+1][nn] = (-sa - w * a[i][nn] - q * a[i][nn-1]) / x;
//                            }
//                            else {
//                                v = cdiv(compl(-r-y*a[i][nn-1], -s-y*a[i][nn]),compl(z,q));
//                                a[i+1][nn-1] = v.re;
//                                a[i+1][nn] = v.im;
//                            }
//                        }
//                    }
//                }
//            }
//            else if (q == 0) {                             /* real vector */
//                m = nn;
//                a[nn][nn] = 1;
//                for (i = nn - 1; i >= 0; i--) {
//                    w = a[i][i] - p;
//                    r = a[i][nn];
//                    for (j = m; j < nn; j++) {
//                        r += a[i][j] * a[j][nn];
//                    }
//                    if (wi[i] < 0) {
//                        z = w;
//                        s = r;
//                    }
//                    else {
//                        m = i;
//                        if (wi[i] == 0) {
//                            if ((t = w) == 0) t = eps * anorm;
//                            a[i][nn] = -r / t;
//                        }
//                        else {            /* solve real equations */
//                            x = a[i][i+1];
//                            y = a[i+1][i];
//                            q = (wr[i] - p) * (wr[i] - p) + wi[i]*wi[i];
//                            t = (x * s - z * r) / q;
//                            a[i][nn] = t;
//                            if (fabs(x) <= fabs(z)) {
//                                a[i+1][nn] = (-s - y * t) / z;
//                            }
//                            else {
//                                a[i+1][nn] = (-r - w * t) / x;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        /* vectors of isolated roots */
//        for (i = 0; i < n; i++) {
//            if (i < low || i > hi-1) {
//                for (j = i; j < n; j++) {
//                    vr[i][j] = a[i][j];
//                }
//            }
//        }
//        /* multiply by transformation matrix */
//        
//        for (j = n-1; j >= low; j--) {
//            m = IMIN(j,hi-1);
//            for (i = low; i <= nn; i++) {
//                for (z = 0,k = low; k <= m; k++) {
//                    z += vr[i][k] * a[k][j];
//                }
//                vr[i][j] = z;
//            }
//        }
//    }
//    /* rearrange complex eigenvectors */
//    for (j = 0; j < n; j++) {
//        if (wi[j] != 0) {
//            for (i = 0; i < n; i++) {
//                vi[i][j] = vr[i][j+1];
//                vr[i][j+1] = vr[i][j];
//                vi[i][j+1] = -vi[i][j];
//            }
//            j++;
//        }
//    }
//}
