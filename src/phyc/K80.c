/*
 *  K80.c
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

#include "K80.h"

#include "matrix.h"

/***************
 K80
 R matrix
 
    A C G T
 A  * 1 k 1
 C  1 * 1 k
 G  k 1 * 1
 T  1 k 1 *
 
 freqs  A C G T = 0.25
 0 1 2 3
 
 rates a b c d e
 0 1 2 3 4
 
 Normalized Q = R.pi
      -1 1/(k+2) k/(k+2) 1/(k+2)
 1/(k+2)     -1  1/(k+2) k/(k+2)
 k/(k+2) 1/(k+2)      -1 1/(k+2)
 1/(k+2) k/(k+2) 1/(k+2)      -1
 
 Eigen values:
 [0, -(2k+2)/(k+2), -(2k+2)/(k+2), -4/(k+2)]
 
 Eigen vector U
 1/4  1/4  1/4  1/4
   0  1/2    0 -1/2
 1/2    0 -1/2    0
 1/4 -1/4  1/4 -1/4

 Eigen vector U^-1
 1  0  1  1
 1  1  0 -1
 1  0 -1  1
 1 -1  0 -1
 
 Derivative of Q with respect to k
 
          0 -1/(k+2)^2  2/(k+2)^2 -1/(k+2)^2
 -1/(k+2)^2          0 -1/(k+2)^2  2/(k+2)^2
  2/(k+2)^2 -1/(k+2)^2          0 -1/(k+2)^2
 -1/(k+2)^2  2/(k+2)^2 -1/(k+2)^2          0
 ***************/

static void k80_update_Q( SubstitutionModel *m );

SubstitutionModel * new_K80(){
	Simplex* freqs = new_Simplex("k80", 4);
	for (int i = 0; i < Parameters_count(freqs->parameters); i++) {
		Parameters_set_estimate(freqs->parameters, false, i);
	}
    SubstitutionModel *m = create_nucleotide_model("K80", K80, freqs);
    
    m->update_Q = k80_update_Q;
    
    m->rates = new_Parameters( 1 );
    Parameters_move(m->rates, new_Parameter_with_postfix("k80.kappa", "model", 3, new_Constraint(0.0001, 100) ) );
    
    return m;
}

SubstitutionModel * new_K80_with_values( const double kappa ){
	Parameter * k = new_Parameter_with_postfix("k80.kappa", "model", kappa, new_Constraint(0.0001, 100));
	SubstitutionModel* m = new_K80_with_parameters(k);
	free_Parameter(k);
	return m;
}

SubstitutionModel * new_K80_with_parameters( Parameter* kappa ){
	Simplex* freqs = new_Simplex("k80", 4);
	for (int i = 0; i < Parameters_count(freqs->parameters); i++) {
		Parameters_set_estimate(freqs->parameters, false, i);
	}
	SubstitutionModel *m = create_nucleotide_model("K80", K80, freqs);
	
	m->update_Q = k80_update_Q;
	
	m->rates = new_Parameters( 1 );
	Parameters_add(m->rates, kappa);
	
	return m;
}

void k80_update_Q( SubstitutionModel *m ){
	
	//    m->Q[1][3] = m->Q[3][1] = m->Q[0][2] = m->Q[2][0] = Parameters_value(m->rates, 0)*0.25; // kappa
	//    m->Q[0][1] = m->Q[1][0] = m->Q[0][3] = m->Q[3][0] = m->Q[1][2] = m->Q[2][1] = m->Q[2][3] = m->Q[3][2] = 0.25;
	//update_eigen_system( m );
	m->Q[1][3] = m->Q[3][1] = m->Q[0][2] = m->Q[2][0] = Parameters_value(m->rates, 0)/(Parameters_value(m->rates, 0)+2.0);
	m->Q[0][1] = m->Q[1][0] = m->Q[0][3] = m->Q[3][0] = m->Q[1][2] = m->Q[2][1] = m->Q[2][3] = m->Q[3][2] = 1.0/(Parameters_value(m->rates, 0)+2.0);
	m->Q[0][0] = m->Q[1][1] = m->Q[2][2] = m->Q[3][3] = -1;
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    m->need_update = false;
}

void k80_dQ(SubstitutionModel *m, int index, double* mat, double t){
	if( m->need_update ){
		m->update_Q(m);
	}
	
	const double k = Parameters_value(m->rates, 0);
	const double iev[16] = {1,0,1,1,1,1,0,-1,1,0,-1,1,1,-1,0,-1};
	const double ev[16] = {0.25,0.25,0.25,0.25,0,0.5,0,-0.5,0.5,0,-0.5,0,0.25,-0.25,0.25,-0.25};
	const double *v = m->eigendcmp->eval;
	//	const double v[4] = {0, -(2*k+2)/(k+2), -(2*k+2)/(k+2), -4/(k+2)};
	double x[16];
	x[0] = x[5] = x[10] = x[15] = 0;
	x[1] = x[3] = x[4] = x[6] = x[9] = x[11] = x[12] = x[14] = -1.0/pow(k+2.0, 2);
	x[2] = x[7] = x[8] = x[13]= 2.0/pow(k+2.0, 2);
	
	//		Matrix_mult2(mat, iev, x, 4,4,4,4);
	//		Matrix_mult2(x, mat, ev, 4,4,4,4);
	Matrix_mult3(mat, (const double**)m->eigendcmp->Invevec, x, 4,4,4,4);
	Matrix_mult4(x, mat, (const double**)m->eigendcmp->evec, 4,4,4,4);
	// up to now the above operations can be recycled across branches
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			if(v[i] != v[j]){
				mat[i*4+j] = x[i*4+j]*(exp(v[i]*t) - exp(v[j]*t))/(v[i]-v[j]);
			}
			else{
				mat[i*4+j] = x[i*4+j]*t*exp(v[i]*t);
			}
		}
	}
	//		Matrix_mult2(x, ev, mat, 4,4,4,4);
	//		Matrix_mult2(mat, x, iev, 4,4,4,4);
	Matrix_mult3(x, (const double**)m->eigendcmp->evec, mat, 4,4,4,4);
	Matrix_mult4(mat, x, (const double**)m->eigendcmp->Invevec, 4,4,4,4);
}

