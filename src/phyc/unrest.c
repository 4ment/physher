/*
 *  unrest.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/12/2015.
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

#include "unrest.h"

#include "matrix.h"

static void _nuc_unrestricted_update_Q( SubstitutionModel *m );

int matexp (double **Q, double *P, double t, int n, int TimeSquare );

// row major
static void _p_t_unrestricted( SubstitutionModel *m, const double t, double *P ){
    
    if( m->need_update ){
        m->update_Q(m);
    }
    
    matexp (m->Q, P, t, m->nstate, 10);
}

// t(P D P-1) = t(P-1) D t(P)
// or column major
static void _p_t_transpose_unrestricted( SubstitutionModel *m, const double t, double *P ){
    
    if( m->need_update ){
        m->update_Q(m);
    }
    int n = m->nstate;
    
    matexp (m->Q, P, t, m->nstate, 10);
    
    double v;
    for ( int i = 0; i < n; i++ ){
        for ( int j = 0; j < i; j++ ){
            if ( i !=j ) {
                v = P[i*n+j];
                P[i*n+j] = P[j*n+i];
                P[j*n+i] = v;
            }
        }
    }
    
}


SubstitutionModel * new_UnrestrictedNucleotideModel( ){
	Simplex* freqs = new_Simplex("urev", 4);
	SubstitutionModel *m = create_nucleotide_model("UREV", NON_REVERSIBLE_DNA, freqs);
	
	m->rates = new_Parameters( 11 );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r1",  "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r2",  "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r3",  "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r4",  "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r5",  "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r6",  "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r7",  "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r8",  "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r9",  "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r10", "model", 1, new_Constraint(0.001, 100) ) );
	Parameters_move(m->rates, new_Parameter_with_postfix("unres.r11", "model", 1, new_Constraint(0.001, 100) ) );
	
	
	m->update_Q = _nuc_unrestricted_update_Q;
	
	//m->p_t = _p_t_unrestricted;
	//m->p_t_transpose = _p_t_transpose_unrestricted;
	
	
	return m;
}

SubstitutionModel * new_UnrestrictedNucleotideModel_with_parameters(Simplex* freqs, const Parameters* rates ){
	SubstitutionModel *m = create_nucleotide_model("UREV", NON_REVERSIBLE_DNA, freqs);
	
	m->rates = new_Parameters(11);
	for(int i = 0; i < Parameters_count(rates); i++){
		Parameters_add(m->rates, Parameters_at(rates, i) );
	}
	
	m->update_Q = _nuc_unrestricted_update_Q;
	
	//m->p_t = _p_t_unrestricted;
	//m->p_t_transpose = _p_t_transpose_unrestricted;
	
	
	return m;
}



int matinv( double x[], int n, int m, double space[] ) {
    /* x[n*m]  ... m>=n
     space[n].  This puts the fabs(|x|) into space[0].  Check and calculate |x|.
     Det may have the wrong sign.  Check and fix.
     */
    int i,j,k;
    int *irow = (int*)space;
    double ee = 1e-100, t,t1,xmax, det=1;
    
    for( i = 0; i < n; i++ ) irow[i] = i;
    
    for( i = 0; i < n; i++ )  {
        xmax = fabs(x[i*m+i]);
        for( j = i+1; j < n; j++ ){
            if( xmax < fabs(x[j*m+i]) ){
                xmax = fabs(x[j*m+i]);
                irow[i]=j;
            }
        }
        det *= x[irow[i]*m+i];
        
        if ( xmax < ee )   {
            printf("\nxmax = %.4e close to zero at %3d!\t\n", xmax,i+1);
            exit(-1);
        }
        if( irow[i] != i ) {
            for(j=0; j < m; j++ ) {
                t = x[i*m+j];
                x[i*m+j] = x[irow[i]*m+j];
                x[irow[i]*m+j] = t;
            }
        }
        t = 1./x[i*m+i];
        for( j = 0; j < n; j++ ){
            if (j == i) continue;
            t1 = t*x[j*m+i];
            for(k = 0; k <m; k++ )  x[j*m+k] -= t1*x[i*m+k];
            x[j*m+i] = -t1;
        }
        for( j = 0; j < m; j++)   x[i*m+j] *= t;
        x[i*m+i] = t;
    }
    
    for ( i = n-1; i >= 0; i-- ) {
        if( irow[i] == i) continue;
        for( j = 0; j < n; j++ )  {
            t = x[j*m+i];
            x[j*m+i] = x[j*m + irow[i]];
            x[j*m + irow[i]] = t;
        }
    }
    space[0] = det;
    return(0);
}


// Q[,4]=c(1,1,1,1)
// solve(t(Q),c(0,0,0,1))
// or solve(Q) and take last row as pi
int QtoPi ( double **Q, double *pi, int n ){
    /* from rate matrix Q[] to pi, the stationary frequencies:
     pi * Q = Q' * pi = 0     pi * 1 = 1
     space[] is of size n*(n+1).
     */
    int i,j;
    //double *T = space;      /* T[n*(n+1)]  */
    double T[20];
    
    for( i = 0; i < n+1; i++ ) T[i] = 1;
    
    for( i = 1; i < n; i++ ) {
        for( j = 0; j < n; j++ ){
            T[i*(n+1)+j] =  Q[j][i];     /* transpose */
        }
        T[i*(n+1)+n] = 0.;
    }
    matinv(T, n, n+1, pi);
    
    for( i = 0; i < n; i++ ){
        pi[i] = T[i*(n+1)+n];
    }
    return (0);
}

int matby ( const double a[], const double b[], double c[], int n, int m, int k)
/* a[n*m], b[m*k], c[n*k]  ......  c = a*b
 */
{
    int i,j,i1;
    double t;
    for( i = 0; i < n; i++ ){
        for(j=0; j<k; j++) {
            for (i1=0,t=0; i1<m; i1++){
                t += a[i*m+i1]*b[i1*k+j];
            }
            c[i*k+j] = t;
        }
    }
    return (0);
}

int matexp (double **Q, double *P, double t, int n, int TimeSquare ) {
    /* This calculates the matrix exponential P(t) = exp(t*Q).
     Input: Q[] has the rate matrix, and t is the time or branch length.
     TimeSquare is the number of times the matrix is squared and should
     be from 5 to 31.
     Output: Q[] has the transition probability matrix, that is P(Qt).
     space[n*n]: required working space.
     
     P(t) = (I + Qt/m + (Qt/m)^2/2)^m, with m = 2^TimeSquare.
     
     See equation (2.22) in Yang (2006) and the discussion below it.
     T[it=0] is the current matrix, and T[it=1] is the squared result matrix,
     used to avoid copying matrices.
     Use an even TimeSquare to avoid one round of matrix copying.
     */
    int it, i, j;
    double *T[2];
    
    if(TimeSquare<2 || TimeSquare>31) error("TimeSquare not good");
    
    T[0] = dvector(n*n);
    T[1] = P;
    
    // T[0] = Qt/m
    for(i = 0; i < n; i++ ){
        for(j = 0; j < n; j++ ){
            T[0][i*n+j] = ldexp( Q[i][j]*t, -TimeSquare );
        }
    }
    
    // T[1] = T[0] * T[0] = Qt/m * Qt/m
    matby (T[0], T[0], T[1], n, n, n);
    
    // T[0] = T[0] + T1/2 = Qt/m + (Qt/m)^2/2
    for( i = 0; i < n*n; i++ ){
        T[0][i] += T[1][i]/2;
    }
    
    // T[0] = I + T[0] + T1/2 = I + Qt/m + (Qt/m)^2/2
    for( i = 0; i < n; i++ ){
        T[0][i*n+i]++;
    }
    
    for( i = 0, it = 0; i < TimeSquare; i++ ) {
        it = !it;
        matby (T[1-it], T[1-it], T[it], n, n, n);
    }
    if( it == 0 ){
        memcpy(P, T[0], n*n*sizeof(double));
    }
    free(T[0]);
    return(0);
}

void _nuc_unrestricted_update_Q( SubstitutionModel *m ){
    int index = 0;
    for ( int i = 0; i < 4; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
            if( i != j ){
                if( index == 11 ){
                    m->Q[i][j] = 1;
                }
                else{
                    m->Q[i][j] = Parameters_value(m->rates, index);
                    index++;
                }
            }
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
	double f[4];
    QtoPi( m->Q, f, m->nstate);
	m->simplex->set_values(m->simplex, f);
    normalize_Q( m->Q, f, m->nstate );
    
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    
    m->need_update = false;
}
