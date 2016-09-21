/*
 *  hessenberg.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/7/10.
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

#include "hessenberg.h"
#include <math.h>
#include <string.h>

#include "utils.h"
#include "matrix.h"


#define SWAP(g,h){y=(g);(g)=(h);(h)=y;}
#define RADIX 2.0
#define TINY 1.0e-20


typedef struct complex{
	double re;
	double im;
}complex;

complex compl (double re,double im){
    complex r;
	
    r.re = re;
    r.im = im;
    return(r);
}


complex cdiv (complex a,complex b){
    double ratio, den;
    complex c;
	
    if (fabs(b.re) <= fabs(b.im)) {
        ratio = b.re / b.im;
        den = b.im * (1 + ratio * ratio);
        c.re = (a.re * ratio + a.im) / den;
        c.im = (a.im * ratio - a.re) / den;
    }
    else {
        ratio = b.im / b.re;
        den = b.re * (1 + ratio * ratio);
        c.re = (a.re + a.im * ratio) / den;
        c.im = (a.im - a.re * ratio) / den;
    }
    return(c);
}

// Helper variables for mcdiv
double cr, ci;

void mcdiv(double ar, double ai, double br, double bi) {
	double s, ars, ais, brs, bis;
	
	s = fabs(br) + fabs(bi);
	ars = ar / s;
	ais = ai / s;
	brs = br / s;
	bis = bi / s;
	s = brs * brs + bis * bis;
	cr = (ars * brs + ais * bis) / s;
	ci = (ais * brs - ars * bis) / s;
}


// Given a matrix a[0..n-1][0..n-1], this routine replaces it by a balanced matrix with identical eigenvalues.
// A symmetric matrix is already balanced and is unaffected by this procedure.
// The parameter RADIX should be the machine’s floating-point radix.
void balance(double **a, int n){
	int last,j,i;
	double s,r,g,f,c,sqrdx;
	sqrdx = RADIX*RADIX;
	last=0;
	
	while(last == 0) {
		last=1;
		//Calculate row and column norms.
		for ( i = 1 ;i < n; i++ ) {
			r=c=0.0;
			for ( j = 1 ; j < n; j++ ){
				if (j != i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			}
			
			// If both are nonzero
			if (c && r){
				g=r/RADIX;
				f=1.0;
				s=c+r;
				// Find the integer power of the machine radix that comes closest to balancing the matrix. 
				while (c<g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g=r*RADIX;
				while (c>g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s) {
					last=0;
					g=1.0/f;
					for (j=1;j<=n;j++) a[i][j] *= g;
					for (j=1;j<=n;j++) a[j][i] *= f;
				}
			}
		}
	}
}

/*
 * Transform back eigenvectors of a balanced matrix
 * into the eigenvectors of the original matrix
 */
void unbalance(int n,double **vr,double **vi, int low, int hi, double *scale){
    int i,j,k;
    double tmp;
	
    for (i = low; i <= hi; i++) {
        for (j = 0; j < n; j++) {
            vr[i][j] *= scale[i];
            vi[i][j] *= scale[i];
        }
    }
	
    for (i = low - 1; i >= 0; i--) {
        if ((k = (int)scale[i]) != i) {
            for (j = 0; j < n; j++) {
                tmp = vr[i][j];
                vr[i][j] = vr[k][j];
                vr[k][j] = tmp;
				
                tmp = vi[i][j];
                vi[i][j] = vi[k][j];
                vi[k][j] = tmp;        
            }
        }
    }
	
    for (i = hi + 1; i < n; i++) {
        if ((k = (int)scale[i]) != i) {
            for (j = 0; j < n; j++) {
                tmp = vr[i][j];
                vr[i][j] = vr[k][j];
                vr[k][j] = tmp;
				
                tmp = vi[i][j];
                vi[i][j] = vi[k][j];
                vi[k][j] = tmp;        
            }
        }
    }
}

void eltran(double **a, double **zz, int *order, const int n) {
	int i, j, m;
	
	for ( i = 0; i < n; i++ ) {
		for ( j = i + 1; j < n; j++ ) {
			zz[i][j] = 0.0;
			zz[j][i] = 0.0;
		}
		zz[i][i] = 1.0;
	}
	if ( n <= 2 ) {
		return;
	}
	for ( m = n - 1; m >= 2; m-- ) {
		for (i = m; i < n; i++) {
			zz[i][m - 1] = a[i][m - 2];
		}
		i = order[m - 1];
		if ( i != m ) {
			for ( j = m - 1; j < n; j++ ) {
				zz[m - 1][j] = zz[i - 1][j];
				zz[i - 1][j] = 0.0;
			}
			zz[i - 1][m - 1] = 1.0;
		}
	}
}

// Reduction to Hessenberg form by the elimination method.
// The real, nonsymmetric matrix a[0..n-1][0..n-1] is replaced by an upper Hessenberg matrix with identical eigenvalues. Recommended, but not required, is that this routine be preceded by balanc.
// On output, the Hessenberg matrix is in elements a[i][j] with i ≤ j+1. Elements with i > j+1 are to be thought of as zero, but are returned with random values.

void elmhes(double **a, int *ordr, int n){
	int m, j, i;
	double y, x;
	
	for (i = 0; i < n; i++) ordr[i] = 0;
	
	// m is called r + 1 in the text.
	for (m = 2; m < n; m++) {
		x = 0.0;
		i = m;
		// Find the pivot.
		for (j = m; j <= n; j++) {
			if (fabs(a[j - 1][m - 2]) > fabs(x)) {
				x = a[j - 1][m - 2];
				i = j;
			}
		}
		//Interchange rows and columns.
		ordr[m - 1] = i;
		if (i != m) {
			for (j = m - 2; j < n; j++) {
				y = a[i - 1][j];
				a[i - 1][j] = a[m - 1][j];
				a[m - 1][j] = y;
			}
			for (j = 0; j < n; j++) {
				y = a[j][i - 1];
				a[j][i - 1] = a[j][m - 1];
				a[j][m - 1] = y;
			}
		}
		// Carry out the elimination.
		if (x != 0.0) {
			for (i = m; i < n; i++)	{
				y = a[i][m - 2];
				if (y != 0.0){
					y /= x;
					a[i][m - 2] = y;
					for (j = m - 1; j < n; j++) {
						a[i][j] -= y * a[m - 1][j];
					}
					for (j = 0; j < n; j++) {
						a[j][m - 1] += y * a[j][i];
					}
				}
			}
		}
	}
}

/**
 Nonsymmetric reduction to Hessenberg form by orthogonal similarity transformations.
 */
void orthes ( int n, double **H, double **V, double *ort) {
    //  This is derived from the Algol procedures orthes and ortran,
    //  by Martin and Wilkinson, Handbook for Auto. Comp.,
    //  Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutines in EISPACK.
    
    int low = 0;
    int high = n-1;
    
    for (int m = low+1; m <= high-1; m++) {
        
        // Scale column.
        
        double scale = 0.0;
        for (int i = m; i <= high; i++) {
            scale = scale + fabs(H[i][m-1]);
        }
        if (scale != 0.0) {
            
            // Compute Householder transformation.
            
            double h = 0.0;
            for (int i = high; i >= m; i--) {
                ort[i] = H[i][m-1]/scale;
                h += ort[i] * ort[i];
            }
            double g = sqrt(h);
            if (ort[m] > 0) {
                g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;
            
            // Apply Householder similarity transformation
            // H = (I-u*u'/h)*H*(I-u*u')/h)
            
            for (int j = m; j < n; j++) {
                double f = 0.0;
                for (int i = high; i >= m; i--) {
                    f += ort[i]*H[i][j];
                }
                f = f/h;
                for (int i = m; i <= high; i++) {
                    H[i][j] -= f*ort[i];
                }
            }
            
            for (int i = 0; i <= high; i++) {
                double f = 0.0;
                for (int j = high; j >= m; j--) {
                    f += ort[j]*H[i][j];
                }
                f = f/h;
                for (int j = m; j <= high; j++) {
                    H[i][j] -= f*ort[j];
                }
            }
            ort[m] = scale*ort[m];
            H[m][m-1] = scale*g;
        }
    }
    
    // Accumulate transformations (Algol's ortran).
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            V[i][j] = (i == j ? 1.0 : 0.0);
        }
    }
    
    for (int m = high-1; m >= low+1; m--) {
        if (H[m][m-1] != 0.0) {
            for (int i = m+1; i <= high; i++) {
                ort[i] = H[i][m-1];
            }
            for (int j = m; j <= high; j++) {
                double g = 0.0;
                for (int i = m; i <= high; i++) {
                    g += ort[i] * V[i][j];
                }
                // Double division avoids possible underflow
                g = (g / ort[m]) / H[m][m-1];
                for (int i = m; i <= high; i++) {
                    V[i][j] += g * ort[i];
                }
            }
        }
    }
}

// Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n].
// On input a can be exactly as output from elmhes §11.5; on output it is destroyed.
// The real and imaginary parts of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
void hqr(double **a, int n, double wr[], double wi[]){
	int nn,m,l,k,j,its,i,mmin;
	double z,y,x,w,v,u,t,s,r,q,p,anorm;
	anorm = p = q = r = 0.0;
	
	// Compute matrix norm for possible use in locating single small subdiagonal element.
	for ( i = 0; i < n; i++ )
		for ( j = imax(i-1,0); j < n; j++ )
			anorm += fabs(a[i][j]);
	
	nn = n-1;
	t=0.0;
	// Gets changed only by an exceptional shift.
	// Begin search for next eigenvalue.
	while(nn >= 0) {
		its=0;
		do{
			// Begin iteration: look for single small subdiagonal element.
			for( l = nn; l >= 1; l-- ) {
				s = fabs(a[l-1][l-1]) + fabs(a[l][l]);
				if (s == 0.0) s = anorm;
				if ( fabs(a[l][l-1] + s) == s) break;
			}
			x = a[nn][nn];
			
			//One root found.
			if (l == nn) {
				wr[nn] = x+t;
				wi[nn--] = 0.0;
			}
			else {
				y=a[nn-1][nn-1];
				w=a[nn][nn-1]*a[nn-1][nn];
				//Two roots found...
				if (l == (nn-1)) {
					p = 0.5*(y-x);
					q = p*p+w;
					z = sqrt(fabs(q));
					x += t;
					//	...a real pair.
					if (q >= 0.0) { 
						z =  p + SIGN(z,p);
						wr[nn-1] = wr[nn] = x+z;
						if (z) wr[nn] = x-w/z;
						wi[nn-1] = wi[nn]=0.0;
					} 
					//...a complex pair.
					else{
						wr[nn-1] = wr[nn]=x+p; 
						wi[nn-1] = -(wi[nn]=z);
					}
					nn -= 2;
				}
				// No roots found. Continue iteration.
				else{
					if (its == 30) error("Too many iterations in hqr");
					// Form exceptional shift.
					if (its == 10 || its == 20) {
						t += x;
						for ( i=0; i<=nn; i++ ) a[i][i] -= x;
						s = fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
						y = x = 0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					//Form shift and then look for 2 consecutive small sub-diagonal elements.
					for( m=(nn-2); m>=l; m-- ) {
						z = a[m][m];
						r = x-z;
						s = y-z;
						p = (r*s-w)/a[m+1][m]+a[m][m+1]; //Equation (11.6.23).
						q = a[m+1][m+1]-z-r-s;
						r = a[m+2][m+1];
						
						
						// Scale to prevent overflow or underflow.
						s=fabs(p)+fabs(q)+fabs(r);
						
						p /= s; 
						q /= s; 
						r /= s; 
						if (m == l) break;
						u = fabs(a[m][m-1])*(fabs(q)+fabs(r));
						v = fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1])); //Equation (11.6.26).
						if( (double)(u+v) == v) break;
					}
					for ( i = m+2; i <= nn; i++ ) {
						a[i][i-2]=0.0; 
						if (i != (m+2)) a[i][i-3]=0.0;
					} 
					for ( k = m; k <= nn-1; k++ ) {
						// Double QR step on rows l to nn and columns m to nn.
						if (k != m) {
							p = a[k][k-1];	//Begin setup of Householder vector
							q = a[k+1][k-1];
							r=0.0;
							if (k != (nn-1)) r=a[k+2][k-1];
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
								//Scale to prevent overflow or underflow.
								p /= x;	
								q /= x;
								r /= x;
							}
						} 
						if ((s = SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							if (k == m) {
								if (l != m) 
									a[k][k-1] = -a[k][k-1];
							}
							else a[k][k-1] = -s*x; 
							p += s;     //Equations (11.6.24).
							x = p/s;
							y = q/s;
							z = r/s;
							q /= p;
							r /= p;
							
							// Row modification.
							for ( j = k; j <= nn; j++ ) {
								p = a[k][j]+q*a[k+1][j];
								if (k != (nn-1)) {
									p += r*a[k+2][j]; 
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j]   -= p*x;
							} 
							mmin = nn < k+3 ? nn : k+3;
							// Column modification.
							for ( i = l; i <= mmin; i++ ) {
								p=x*a[i][k]+y*a[i][k+1];
								if (k != (nn-1)) {
									p += z*a[i][k+2]; 
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q; 
								a[i][k]   -= p;
							}
						}
					}
				}
			}
		}
		while( 1 < nn-1);
	}
}




void hqr3(int n, int low, int hgh, double **h, double **zz, double *wr, double *wi) {
	int i, j, k, l = 0, m, en, na, itn, its;
	double p = 0, q = 0, r = 0, s = 0, t, w, x = 0, y, ra, sa, vi, vr, z = 0, norm, tst1, tst2;
	bool notLast;
	
	
	norm = 0.0;
	k = 1;
	/* store isolated roots and compute matrix norm */
	for (i = 0; i < n; i++) {
		for (j = k - 1; j < n; j++) {
			norm += fabs(h[i][j]);
		}
		k = i + 1;
		if (i + 1 < low || i + 1 > hgh) {
			wr[i] = h[i][i];
			wi[i] = 0.0;
		}
	}
	en = hgh;
	t = 0.0;
	itn = n * 30;
	while (en >= low) {    /* search for next eigenvalues */
		its = 0;
		na = en - 1;
		while (en >= 1) {
			/* look for single small sub-diagonal element */
			bool fullLoop = true;
			for (l = en; l > low; l--) {
				s = fabs(h[l - 2][l - 2]) + fabs(h[l - 1][l - 1]);
				if (s == 0.0) {
					s = norm;
				}
				tst1 = s;
				tst2 = tst1 + fabs(h[l - 1][l - 2]);
				if (tst2 == tst1) {
					fullLoop = false;
					break;
				}
			}
			if (fullLoop) {
				l = low;
			}
			
			x = h[en - 1][en - 1];    /* form shift */
			if (l == en || l == na) {
				break;
			}
			if (itn == 0) {
				/* eigenvalues have not converged */
				error("Eigenvalues not converged");
			}
			y = h[na - 1][na - 1];
			w = h[en - 1][na - 1] * h[na - 1][en - 1];
			/* form exceptional shift */
			if (its == 10 || its == 20) {
				t += x;
				for (i = low - 1; i < en; i++) {
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
			for (m = en - 2; m >= l; m--) {
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
				if (m == l) {
					break;
				}
				tst1 = fabs(p) * (fabs(h[m - 2][m - 2]) + fabs(z) + fabs(h[m][m]));
				tst2 = tst1 + fabs(h[m - 1][m - 2]) * (fabs(q) + fabs(r));
				if (tst2 == tst1) {
					break;
				}
			}
			for (i = m + 2; i <= en; i++) {
				h[i - 1][i - 3] = 0.0;
				if (i != m + 2) {
					h[i - 1][i - 4] = 0.0;
				}
			}
			for (k = m; k <= na; k++) {
				notLast = k != na;
				if (k != m) {
					p = h[k - 1][k - 2];
					q = h[k][k - 2];
					r = 0.0;
					if (notLast) {
						r = h[k + 1][k - 2];
					}
					x = fabs(p) + fabs(q) + fabs(r);
					if (x != 0.0) {
						p /= x;
						q /= x;
						r /= x;
					}
				}
				if (x != 0.0) {
					if (p < 0.0) {    /* sign */
						s = -sqrt(p * p + q * q + r * r);
					} else {
						s = sqrt(p * p + q * q + r * r);
					}
					if (k != m) {
						h[k - 1][k - 2] = -s * x;
					} else if (l != m) {
						h[k - 1][k - 2] = -h[k - 1][k - 2];
					}
					p += s;
					x = p / s;
					y = q / s;
					z = r / s;
					q /= p;
					r /= p;
					if (!notLast) {
						for (j = k - 1; j < n; j++) {    /* row modification */
							p = h[k - 1][j] + q * h[k][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for (i = 0; i < j; i++) {    /* column modification */
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
					} else {
						for (j = k - 1; j < n; j++) {    /* row modification */
							p = h[k - 1][j] + q * h[k][j] + r * h[k + 1][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
							h[k + 1][j] -= p * z;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for (i = 0; i < j; i++) {    /* column modification */
							p = x * h[i][k - 1] + y * h[i][k] + z * h[i][k + 1];
							h[i][k - 1] -= p;
							h[i][k] -= p * q;
							h[i][k + 1] -= p * r;
						}
						/* accumulate transformations */
						for (i = low - 1; i < hgh; i++) {
							p = x * zz[i][k - 1] + y * zz[i][k] +
							z * zz[i][k + 1];
							zz[i][k - 1] -= p;
							zz[i][k] -= p * q;
							zz[i][k + 1] -= p * r;
						}
					}
				}
			}                 /* for k */
		}                     /* while infinite loop */
		if (l == en) {                 /* one root found */
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
		if (q >= 0.0) {     /* real pair */
			if (p < 0.0) {    /* sign */
				z = p - fabs(z);
			} else {
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
			for (j = na - 1; j < n; j++) {    /* row modification */
				z = h[na - 1][j];
				h[na - 1][j] = q * z + p * h[en - 1][j];
				h[en - 1][j] = q * h[en - 1][j] - p * z;
			}
			for (i = 0; i < en; i++) {    /* column modification */
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
		} else {    /* complex pair */
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
						} else {
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
							} else {    /* solve real equations */
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
			} else if (q > 0.0) {
				m = na;
				if (fabs(h[en - 1][na - 1]) > fabs(h[na - 1][en - 1])) {
					h[na - 1][na - 1] = q / h[en - 1][na - 1];
					h[na - 1][en - 1] = (p - h[en - 1][en - 1]) / h[en - 1][na - 1];
				} else {
					mcdiv(0.0, -h[na - 1][en - 1], h[na - 1][na - 1] - p, q);
					h[na - 1][na - 1] = cr;
					h[na - 1][en - 1] = ci;
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
						} else {
							m = i + 1;
							if (wi[i] == 0.0) {
								mcdiv(-ra, -sa, w, q);
								h[i][na - 1] = cr;
								h[i][en - 1] = ci;
							} else {    /* solve complex equations */
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
								mcdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
								h[i][na - 1] = cr;
								h[i][en - 1] = ci;
								if (fabs(x) > fabs(z) + fabs(q)) {
									h[i + 1]
									[na - 1] = (q * h[i][en - 1] -
                                                w * h[i][na - 1] - ra) / x;
									h[i + 1][en - 1] = (-sa - w * h[i][en - 1] -
														q * h[i][na - 1]) / x;
								} else {
									mcdiv(-r - y * h[i][na - 1], -s - y * h[i][en - 1], z, q);
									h[i + 1][na - 1] = cr;
									h[i + 1][en - 1] = ci;
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
		for (i = 0; i < n; i++) {
			if (i + 1 < low || i + 1 > hgh) {
				for (j = i; j < n; j++) {
					zz[i][j] = h[i][j];
				}
			}
		}
		/* multiply by transformation matrix to give vectors of
		 * original full matrix. */
		for (j = n - 1; j >= low - 1; j--) {
			m = ((j + 1) < hgh) ? (j + 1) : hgh; /* min */
			for (i = low - 1; i < hgh; i++) {
				z = 0.0;
				for (k = low - 1; k < m; k++) {
					z += zz[i][k] * h[k][j];
				}
				zz[i][j] = z;
			}
		}
	}
}

// Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU decomposition of a rowwise permutation of itself.
// a and n are input.
// a is output, arranged as in equation (2.3.14) above;
// indx[0..n-1] is an output vector that records the row permutation effected by the partial pivoting;
// d is output as ±1 depending on whether the number of row interchanges was even or odd, respectively.
// This routine is used in combination with lubksb to solve linear equations or invert a matrix.
void ludcmp(double **a, int n, int *indx, double *d){
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;
	imax = 0;
    
	vv = dvector(n); //vv stores the implicit scaling of each row.
	*d = 1.0;
	
	for( i = 0; i < n; i++ ) {
		big = 0.0;
		for( j = 0; j < n; j++ )
			// No row interchanges yet. Loop over rows to get the implicit scaling information.
			if ((temp=fabs(a[i][j])) > big) big = temp;
		if (big == 0.0) error("Singular matrix in routine ludcmp"); //No nonzero largest element.
		vv[i] = 1.0/big;
		// Save the scaling. This is the loop over columns of Crout’s method.
	}
	for ( j = 0; j < n ;j++ ) {
		// This is equation (2.3.12) except for i = j.
		for ( i = 0; i < j; i++ ) {
			sum = a[i][j];
			for ( k = 0; k < i; k++ ) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for ( i = j; i < n; i++ ) {
			sum=a[i][j];
			for( k = 0; k < j; k++ )
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			// Is the figure of merit for the pivot better than the best so far?
			if( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if(j != imax) { //Do we need to interchange rows?
			for( k = 0; k < n; k++ ) { // Yes, do so...
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax]=vv[j]; //...and change the parity of d. Also interchange the scale factor.
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = TINY;
		//If the pivot element is zero the matrix is singular (at least to the precision of the algorithm). For some applications on singular matrices, it is desirable to substitute TINY for zero.
		if (j != n-1) {	//Now, finally, divide by the pivot element.
			dum = 1.0/(a[j][j]);
			for( i = j+1; i < n; i++ ) a[i][j] *= dum;
			
		}
	}	//Go back for the next column in the reduction.
	
	free(vv);
}

void ludcmp2(double *a, int n, int *indx, double *d){
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;
	imax = 0;
    
	vv = dvector(n); //vv stores the implicit scaling of each row.
	*d = 1.0;
	
	for( i = 0; i < n; i++ ) {
		big = 0.0;
		for( j = 0; j < n; j++ )
			// No row interchanges yet. Loop over rows to get the implicit scaling information.
			if ((temp=fabs(a[i*n+j])) > big) big = temp;
		if (big == 0.0) error("Singular matrix in routine ludcmp"); //No nonzero largest element.
		vv[i] = 1.0/big;
		// Save the scaling. This is the loop over columns of Crout’s method.
	}
	for ( j = 0; j < n ;j++ ) {
		// This is equation (2.3.12) except for i = j.
		for ( i = 0; i < j; i++ ) {
			sum = a[i*n+j];
			for ( k = 0; k < i; k++ ) sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
		}
		big=0.0;
		for ( i = j; i < n; i++ ) {
			sum=a[i*n+j];
			for( k = 0; k < j; k++ )
				sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
			// Is the figure of merit for the pivot better than the best so far?
			if( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if(j != imax) { //Do we need to interchange rows?
			for( k = 0; k < n; k++ ) { // Yes, do so...
				dum = a[imax*n+k];
				a[imax*n+k] = a[j*n+k];
				a[j*n+k] = dum;
			}
			*d = -(*d);
			vv[imax]=vv[j]; //...and change the parity of d. Also interchange the scale factor.
		}
		indx[j] = imax;
		if (a[j*n+j] == 0.0) a[j*n+j] = TINY;
		//If the pivot element is zero the matrix is singular (at least to the precision of the algorithm). For some applications on singular matrices, it is desirable to substitute TINY for zero.
		if (j != n-1) {	//Now, finally, divide by the pivot element.
			dum = 1.0/(a[j*n+j]);
			for( i = j+1; i < n; i++ ) a[i*n+j] *= dum;
			
		}
	}	//Go back for the next column in the reduction.
	
	free(vv);
}

double LUDecompose_det( double **m, const int dim ) {
	double big = 0.0;
	double temp = 0.0;
	double *vv = dvector(dim);
	double d = 1.0;
	
	for ( size_t i = 0; i < dim; i++ ) {
		big = 0.0;
		for ( size_t j = 0; j < dim; j++ ) {
			if( (temp = fabs( m[i][j] )) > big) big = temp;
		}
		
		if ( big == 0.0 ) {
			error("LU decomposition does not like singular matrices");
		}
		vv[i] = 1.0/big;
	}
	
	size_t i, j;
	size_t imax;
	for ( size_t k = 0; k < dim; k++) {
		big = 0.0;
		imax = k;
		for ( i = k+1; i < dim; i++) {
			temp = vv[i] * fabs( m[i][k] );
			if( temp > big ){
				big = temp;
				imax = i;
			}
		}
		
		if( k != imax ){
			for ( j = 0; j < dim; j++) {
				temp = m[imax][j];
				m[imax][j] = m[k][j];
				m[k][j] = temp;
			}
			d = -d;
			vv[imax] = vv[k];
		}
		
		if( m[k][k] == 0.0 ) m[k][k] = 1.0e-25;
		
		for ( i = k+1; i < dim; i++) {
			temp = m[i][k] /= m[k][k];
			for ( j = k+1; j < dim; j++) {
				m[i][j] -= temp * m[k][j];
			}
		}
	}
	
	for ( size_t i = 0; i < dim; i++) {
		d *= m[i][i];
	}
	
	free(vv);
	
	return d;
}

double LUDecompose_det2( double *m, const int dim ) {
	double big = 0.0;
	double temp = 0.0;
	double *vv = dvector(dim);
	double d = 1.0;
	
	for ( size_t i = 0; i < dim; i++ ) {
		big = 0.0;
		for ( size_t j = 0; j < dim; j++ ) {
			if( (temp = fabs( m[i*dim+j] )) > big) big = temp;
		}
		
		if ( big == 0.0 ) {
			error("LU decomposition does not like singular matrices");
		}
		vv[i] = 1.0/big;
	}
	
	size_t i, j;
	size_t imax;
	for ( size_t k = 0; k < dim; k++) {
		big = 0.0;
		imax = k;
		for ( i = k+1; i < dim; i++) {
			temp = vv[i] * fabs( m[i*dim+k] );
			if( temp > big ){
				big = temp;
				imax = i;
			}
		}
		
		if( k != imax ){
			for ( j = 0; j < dim; j++) {
				temp = m[imax*dim+j];
				m[imax*dim+j] = m[k*dim+j];
				m[k*dim+j] = temp;
			}
			d = -d;
			vv[imax] = vv[k];
		}
		
		if( m[k*dim+k] == 0.0 ) m[k*dim+k] = 1.0e-25;
		
		for ( i = k+1; i < dim; i++) {
			temp = m[i*dim+k] /= m[k*dim+k];
			for ( j = k+1; j < dim; j++) {
				m[i*dim+j] -= temp * m[k*dim+j];
			}
		}
	}
	
	for ( size_t i = 0; i < dim; i++) {
		d *= m[i*dim+i];
	}
	free(vv);
	return d;
}

// compute the determinant of m and the compute its inverse
double LUDecompose_det_and_inverse( double *m, const size_t dim ){
	double d = 1.0;
	
	double *y = dvector(dim*dim);
	int *indx = ivector(dim);
	
	//LU decomposition
	double big = 0.0;
	double temp = 0.0;
	double *vv = dvector(dim);
	
	size_t i, j;
	for ( i = 0; i < dim; i++ ) {
		big = 0.0;
		for ( size_t j = 0; j < dim; j++ ) {
			if( (temp = fabs( m[i*dim+j] )) > big) big = temp;
		}
		
		if ( big == 0.0 ) {
			error("LU decomposition does not like singular matrices");
		}
		vv[i] = 1.0/big;
	}
	
	size_t imax;
	for ( size_t k = 0; k < dim; k++) {
		big = 0.0;
		imax = k;
		for ( i = k+1; i < dim; i++) {
			temp = vv[i] * fabs( m[i*dim+k] );
			if( temp > big ){
				big = temp;
				imax = i;
			}
		}
		
		if( k != imax ){
			for ( j = 0; j < dim; j++) {
				temp = m[imax*dim+j];
				m[imax*dim+j] = m[k*dim+j];
				m[k*dim+j] = temp;
			}
			d = -d;
			vv[imax] = vv[k];
		}
		
		indx[k] = imax;
		if( m[k*dim+k] == 0.0 ) m[k*dim+k] = 1.0e-25;
		
		for ( i = k+1; i < dim; i++) {
			temp = m[i*dim+k] /= m[k*dim+k];
			for ( j = k+1; j < dim; j++) {
				m[i*dim+j] -= temp * m[k*dim+j];
			}
		}
	}
	
//	for ( size_t i = 0; i < dim; i++) {
//		d *= m[i*dim+i];
//	}
	
	// determinant
	d = 0;
	for ( i = 0; i < dim; i++) {
		//d *= m[i*dim+i];
		d += log(m[i*dim+i]);
		fprintf(stderr, "LUDecompose_det_and_inverse det: %f\n", d);
	}
	d = exp(d);
	
	// Inverse of matrix
	for ( j = 0 ; j < dim; j++ ) {
		memset( vv, 0.0, dim*sizeof(double));
		vv[j] = 1.0;
		lubksb2(m,dim,indx,vv);
		for ( i = 0; i < dim; i++ ) y[i*dim+j] = vv[i];
    }
	
	memcpy( m, y, dim*dim*sizeof(double) );
	
	free(y);
	free(vv);
	free(indx);
	return d;
}



/* Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix A
 * but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input as 
 * the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector B, 
 * and returns with the solution vector X. a, n, and indx are not modified by this routine and can 
 * be left in place for successive calls with different right-hand sides b. This routine takes into
 * account the possibility that b will begin with many zero elements, so it is efficient for use
 * in matrix inversion. */
void lubksb( double **a, int n, int *indx, double *b ){
	int i,ii=-1,ip,j; 
	double sum;
	//When ii is set to a positive value, it will become the index of the first nonvanishing element of b. We now do the forward substitution, equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go.
	for ( i = 0; i < n; i++ ) { 
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i]; 
		if (ii >= 0 )
			for( j = ii; j < i; j++ ) sum -= a[i][j]*b[j];
		else if (sum) ii=i; // A nonzero element was encountered, so from now on we will have to do the sums in the loop above.
		b[i] = sum;
	}
	for( i = n-1; i >= 0; i-- ) {
		sum = b[i];  
		for( j = i+1; j < n; j++ ) sum -= a[i][j]*b[j];
		b[i] = sum/a[i][i]; //Store a component of the solution vector X.
	}
	
}

void lubksb2( double *a, int n, int *indx, double *b ){
	int i,ii=-1,ip,j; 
	double sum;
	//When ii is set to a positive value, it will become the index of the first nonvanishing element of b. We now do the forward substitution, equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go.
	for ( i = 0; i < n; i++ ) { 
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i]; 
		if (ii >= 0 )
			for( j = ii; j < i; j++ ) sum -= a[i*n+j]*b[j];
		else if (sum) ii=i; // A nonzero element was encountered, so from now on we will have to do the sums in the loop above.
		b[i] = sum;
	}
	for( i = n-1; i >= 0; i-- ) {
		sum = b[i];  
		for( j = i+1; j < n; j++ ) sum -= a[i*n+j]*b[j];
		b[i] = sum/a[i*n+i]; //Store a component of the solution vector X.
	}
	
}


