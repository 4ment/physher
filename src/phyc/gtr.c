/*
 *  gtr.c
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

#include "gtr.h"

#include "nucsubst.h"
#include "matrix.h"

/***************
 Q matrix
 
 A C G T
 A  * a b c
 C  a * d e
 G  b d * 1
 T  c e 1 *
 
 freqs  A C G T
 0 1 2 3
 
 rates a b c d e
 0 1 2 3 4
  
 a = aa;
 b = (1 - aa) *      bb;
 c = (1 - aa) * (1 - bb) *      cc;
 d = (1 - aa) * (1 - bb) * (1 - cc);
 
 aa = a;
 bb = b /  (1 - aa);
 cc = c / ((1 - aa)*(1 - bb));
 
 ***************/


/*! Initialize a general time reversible matrix.
 *  Return a Q matrix (pi+freq)
 */


SubstitutionModel * new_GTR(Simplex* freqs){
	Parameter* ac = new_Parameter_with_postfix("gtr.a", "model", 1, new_Constraint(0.001, 100) );
	Parameter* ag = new_Parameter_with_postfix("gtr.b", "model", 1, new_Constraint(0.001, 100) );
	Parameter* at = new_Parameter_with_postfix("gtr.c", "model", 1, new_Constraint(0.001, 100) );
	Parameter* cg = new_Parameter_with_postfix("gtr.d", "model", 1, new_Constraint(0.001, 100) );
	Parameter* ct = new_Parameter_with_postfix("gtr.e", "model", 1, new_Constraint(0.001, 100) );
	
	SubstitutionModel* model = new_GTR_with_parameters(freqs, ac, ag, at, cg, ct, NULL);
	
	free_Parameter(ac);
	free_Parameter(ag);
	free_Parameter(at);
	free_Parameter(cg);
	free_Parameter(ct);
	
	return model;
}

SubstitutionModel * new_GTR_with_values( const double *freqs, const double *rates ){
    check_frequencies( freqs, 4 );
    
	Parameter* ac = new_Parameter_with_postfix("gtr.a", "model", 1, new_Constraint(0.001, 100) );
	Parameter* ag = new_Parameter_with_postfix("gtr.b", "model", 1, new_Constraint(0.001, 100) );
	Parameter* at = new_Parameter_with_postfix("gtr.c", "model", 1, new_Constraint(0.001, 100) );
	Parameter* cg = new_Parameter_with_postfix("gtr.d", "model", 1, new_Constraint(0.001, 100) );
	Parameter* ct = new_Parameter_with_postfix("gtr.e", "model", 1, new_Constraint(0.001, 100) );
	
    Simplex* sfreqs = new_Simplex("gtr", 4);
	sfreqs->set_values(sfreqs, freqs);
	
	SubstitutionModel* model = new_GTR_with_parameters(sfreqs, ac, ag, at, cg, ct, NULL);

	free_Parameter(ac);
	free_Parameter(ag);
	free_Parameter(at);
	free_Parameter(cg);
	free_Parameter(ct);
    return model;
}

SubstitutionModel * new_GTR_with_parameters( Simplex* freqs, Parameter* ac, Parameter* ag, Parameter* at, Parameter* cg, Parameter* ct, Parameter* gt){
	
	SubstitutionModel *m = create_nucleotide_model("GTR", GTR, freqs);
	
	m->update_Q = gtr_update_Q;
	
	m->dPdp = gtr_dQ;
	m->dQ = dvector(16);
	m->relativeTo = -1;
	
	m->rates = new_Parameters(5);
	if (ac != NULL) {
		Parameters_add(m->rates, ac);
	}
	else{
		m->relativeTo = 0;
	}
	if (ag != NULL) {
		Parameters_add(m->rates, ag);
	}
	else{
		m->relativeTo = 1;
	}
	if (at != NULL) {
		Parameters_add(m->rates, at);
	}
	else{
		m->relativeTo = 2;
	}
	if (cg != NULL) {
		Parameters_add(m->rates, cg);
	}
	else{
		m->relativeTo = 3;
	}
	if (ct != NULL) {
		Parameters_add(m->rates, ct);
	}
	else{
		m->relativeTo = 4;
	}
	if (gt != NULL) {
		Parameters_add(m->rates, gt);
	}
	else{
		m->relativeTo = 5;
	}
	return m;
}

SubstitutionModel * new_GTR_with_simplexes( Simplex* freqs, Simplex* rates){
	
	SubstitutionModel *m = create_nucleotide_model("GTR", GTR, freqs);
	
	m->update_Q = gtr_simplexes_update_Q;
	
	m->dPdp = gtr_dQ;
	m->dQ = dvector(16);
	
	m->rates_simplex = rates;
	
	return m;
}

void gtr_update_Q( SubstitutionModel *m ){
	double temp;
	int index = 0;
	const double* freqs = m->get_frequencies(m);
	for ( int i = 0; i < 4; i++ )  {
		for ( int j = i + 1; j < 4; j++ ) {
			if(m->relativeTo != index){
				temp = Parameters_value(m->rates, index++);
			}
			else{
				temp = 1;
			}
			m->Q[i][j] = temp * freqs[j];
			m->Q[j][i] = temp * freqs[i];
		}
	}
	update_eigen_system( m );
	m->need_update = false;
}

void gtr_simplexes_update_Q( SubstitutionModel *m ){
	int index = 0;
	const double* freqs = m->simplex->get_values(m->simplex);
	const double* rates = m->rates_simplex->get_values(m->rates_simplex);
	for ( int i = 0; i < 4; i++ )  {
		for ( int j = i + 1; j < 4; j++ ) {
			m->Q[i][j] = rates[index] * freqs[j];
			m->Q[j][i] = rates[index] * freqs[i];
			index++;
		}
	}
	update_eigen_system( m );
	m->need_update = false;
}

void nuc_sym_update_Q_unnormalized( SubstitutionModel *m, double* mat ){
	double temp;
	int index = 0;
	const double* freqs = m->get_frequencies(m);
	for ( int i = 0; i < 4; i++ )  {
		for ( int j = i + 1; j < 4; j++ ) {
			if(m->relativeTo != index){
				temp = Parameters_value(m->rates, index++);
			}
			else{
				temp = 1;
			}
			mat[i*4+j] = temp * freqs[j];
			mat[j*4+i] = temp * freqs[i];
		}
	}
	for (int i = 0; i < 4; i++) {
		double sum = 0.0;
		for (int j = 0; j < 4; j++) {
			if (i != j) {
				sum += mat[i*4+j];
			}
		}
		mat[i*4+i] = -sum;
	}
}

void gtr_dQ(SubstitutionModel *m, int index, double* mat, double t){
	
	double *x = m->dQ;
	
	if(m->need_update){
		m->update_Q(m);
	}
	double* Q = dvector(16);
	
	nuc_sym_update_Q_unnormalized(m, Q);
	
	if(m->dQ_need_update){
		const double a = Parameters_value(m->rates, 0);
		const double b = Parameters_value(m->rates, 1);
		const double c = Parameters_value(m->rates, 2);
		const double d = Parameters_value(m->rates, 3);
		const double e = Parameters_value(m->rates, 4);
		
		const double pi_a = m->simplex->get_value(m->simplex, 0);
		const double pi_c = m->simplex->get_value(m->simplex, 1);
		const double pi_g = m->simplex->get_value(m->simplex, 2);
		const double pi_t = m->simplex->get_value(m->simplex, 3);
		
		const double phi1 = Parameters_value(m->simplex->parameters, 0);
		const double phi2 = Parameters_value(m->simplex->parameters, 1);
		const double phi3 = Parameters_value(m->simplex->parameters, 2);
		
		const double norm = (a*pi_c + b*pi_g + c*pi_t)*pi_a + (a*pi_a + d*pi_g + e*pi_t)*pi_c + (b*pi_a + d*pi_c + pi_t)*pi_g + (c*pi_a + e*pi_c + pi_g)*pi_t;
		const double norm2 = norm*norm;
		
		memset(x, 0, sizeof(double)*16);
		// derivative of norm with respect to a
		if(index == 0){
			x[0] = 2.0*(a*pi_c + b*pi_g + c*pi_t)*pi_a*pi_c/norm2 - pi_c/norm;
			x[1] = -(2.0*a*pi_a*pi_c*pi_c)/norm2 + pi_c/norm;
			x[2] = -(2.0*b*pi_a*pi_c*pi_g)/norm2;
			x[3] = -(2.0*c*pi_a*pi_c*pi_t)/norm2;
			
			x[4] = -2.0*a*pi_a*pi_a*pi_c/norm2 + pi_a/norm;
			x[5] = 2.0*(a*pi_a + d*pi_g + e*pi_t)*pi_a*pi_c/norm2 - pi_a/norm;
			x[6] = -2.0*d*pi_a*pi_c*pi_g/norm2;
			x[7] = -2.0*e*pi_a*pi_c*pi_t/norm2;
			
			x[8] = -2.0*b*pi_a*pi_a*pi_c/norm2;
			x[9] = -2.0*d*pi_a*pi_c*pi_c/norm2;
			x[10] = 2.0*(b*pi_a + d*pi_c + pi_t)*pi_a*pi_c/norm2;
			x[11] = -2.0*pi_a*pi_c*pi_t/norm2;
			
			x[12] = -2.0*c*pi_a*pi_a*pi_c/norm2;
			x[13] = -2.0*e*pi_a*pi_c*pi_c/norm2;
			x[14] = -2.0*pi_a*pi_c*pi_g/norm2;
			x[15] = 2.0*(c*pi_a + e*pi_c + pi_g)*pi_a*pi_c/norm2;
		}
		// derivative of norm with respect to b
		else if(index == 1){
			x[0] = 2.0*(a*pi_c + b*pi_g + c*pi_t)*pi_a*pi_g/norm2 - pi_g/norm;
			x[1] = -2.0*a*pi_a*pi_c*pi_g/norm2;
			x[2] = -2.0*b*pi_a*pi_g*pi_g/norm2 + pi_g/norm;
			x[3] = -2*c*pi_a*pi_g*pi_t/norm2;
			
			x[4] = -2.0*a*pi_a*pi_a*pi_g/norm2;
			x[5] = 2.0*(a*pi_a + d*pi_g + e*pi_t)*pi_a*pi_g/norm2;
			x[6] = -2.0*d*pi_a*pi_g*pi_g/norm2;
			x[7] = -2.0*e*pi_a*pi_g*pi_t/norm2;
			
			x[8] = -2.0*b*pi_a*pi_a*pi_g/norm2 + pi_a/norm;
			x[9] = -2.0*d*pi_a*pi_c*pi_g/norm2;
			x[10] = 2.0*(b*pi_a + d*pi_c + pi_t)*pi_a*pi_g/norm2 - pi_a/norm;
			x[11] = -2.0*pi_a*pi_g*pi_t/norm2;
			
			x[12] = -2*c*pi_a*pi_a*pi_g/norm2;
			x[13] = -2*e*pi_a*pi_c*pi_g/norm2;
			x[14] = -2*pi_a*pi_g*pi_g/norm2;
			x[15] = 2*(c*pi_a + e*pi_c + pi_g)*pi_a*pi_g/norm2;
		}
		// derivative of norm with respect to c
		else if(index == 2){
			x[0] = 2.0*(a*pi_c + b*pi_g + c*pi_t)*pi_a*pi_t/norm2 - pi_t/norm;
			x[1] = -2.0*a*pi_a*pi_c*pi_t/norm2;
			x[2] = -2.0*b*pi_a*pi_g*pi_t/norm2;
			x[3] = -2.0*c*pi_a*pi_t*pi_t/norm2 + pi_t/norm;
			
			x[4] = -2.0*a*pi_a*pi_a*pi_t/norm2;
			x[5] = 2.0*(a*pi_a + d*pi_g + e*pi_t)*pi_a*pi_t/norm2;
			x[6] = -2.0*d*pi_a*pi_g*pi_t/norm2;
			x[7] = -2.0*e*pi_a*pi_t*pi_t/norm2;
			
			x[8] = -2.0*b*pi_a*pi_a*pi_t/norm2;
			x[9] = -2.0*d*pi_a*pi_c*pi_t/norm2;
			x[10] = 2.0*(b*pi_a + d*pi_c + pi_t)*pi_a*pi_t/norm2;
			x[11] = -2.0*pi_a*pi_t*pi_t/norm2;
			
			x[12] = -2.0*c*pi_a*pi_a*pi_t/norm2 + pi_a/norm;
			x[13] = -2.0*e*pi_a*pi_c*pi_t/norm2;
			x[14] = -2.0*pi_a*pi_g*pi_t/norm2;
			x[15] = 2.0*(c*pi_a + e*pi_c + pi_g)*pi_a*pi_t/norm2 - pi_a/norm;
		}
		// derivative of norm with respect to d
		else if(index == 3){
			x[0] = 2.0*(a*pi_c + b*pi_g + c*pi_t)*pi_c*pi_g/norm2;
			x[1] = -2.0*a*pi_c*pi_c*pi_g/norm2;
			x[2] = -2.0*b*pi_c*pi_g*pi_g/norm2;
			x[3] = -2.0*c*pi_c*pi_g*pi_t/norm2;
			
			x[4] = -2.0*a*pi_a*pi_c*pi_g/norm2;
			x[5] = 2.0*(a*pi_a + d*pi_g + e*pi_t)*pi_c*pi_g/norm2 - pi_g/norm;
			x[6] = -2.0*d*pi_c*pi_g*pi_g/norm2 + pi_g/norm;
			x[7] = -2.0*e*pi_c*pi_g*pi_t/norm2;
			
			x[8] = -2.0*b*pi_a*pi_c*pi_g/norm2;
			x[9] = -2.0*d*pi_c*pi_c*pi_g/norm2 + pi_c/norm;
			x[10] = 2.0*(b*pi_a + d*pi_c + pi_t)*pi_c*pi_g/norm2 - pi_c/norm;
			x[11] = -2.0*pi_c*pi_g*pi_t/norm2;
			
			x[12] = -2.0*c*pi_a*pi_c*pi_g/norm2;
			x[13] = -2.0*e*pi_c*pi_c*pi_g/norm2;
			x[14] = -2.0*pi_c*pi_g*pi_g/norm2;
			x[15] = 2.0*(c*pi_a + e*pi_c + pi_g)*pi_c*pi_g/norm2;
		}
		// derivative of norm with respect to e
		else if(index == 4){
			x[0] = 2*(a*pi_c + b*pi_g + c*pi_t)*pi_c*pi_t/norm2;
			x[1] = -2*a*pi_c*pi_c*pi_t/norm2;
			x[2] = -2*b*pi_c*pi_g*pi_t/norm2;
			x[3] = -2*c*pi_c*pi_t*pi_t/norm2;
			
			x[4] = -2*a*pi_a*pi_c*pi_t/norm2;
			x[5] = 2*(a*pi_a + d*pi_g + e*pi_t)*pi_c*pi_t/norm2 - pi_t/norm;
			x[6] = -2*d*pi_c*pi_g*pi_t/norm2;
			x[7] = -2*e*pi_c*pi_t*pi_t/norm2 + pi_t/norm;
			
			x[8] = -2*b*pi_a*pi_c*pi_t/norm2;
			x[9] = -2*d*pi_c*pi_c*pi_t/norm2;
			x[10] = 2*(b*pi_a + d*pi_c + pi_t)*pi_c*pi_t/norm2;
			x[11] = -2*pi_c*pi_t*pi_t/norm2;
			
			x[12] = -2*c*pi_a*pi_c*pi_t/norm2;
			x[13] = -2*e*pi_c*pi_c*pi_t/norm2 + pi_c/norm;
			x[14] = -2*pi_c*pi_g*pi_t/norm2;
			x[15] = 2*(c*pi_a + e*pi_c + pi_g)*pi_c*pi_t/norm2 - pi_c/norm;
		}
		// dQdphi1
		else if(index == 5){
			double dnorm = -4.0*(e*phi1 - e)*phi2*phi2 - 4.0*((phi1 - 1.0)*phi2*phi2 - 2.0*(phi1 - 1.0)*phi2 + phi1 - 1.0)*phi3*phi3 - 4.0*c*phi1 - 2.0*(2.0*(a - c - e)*phi1 - a + c + 2.0*e)*phi2 - 2.0*(2.0*((d - e - 1.0)*phi1 - d + e + 1.0)*phi2*phi2 + 2.0*(b - c - 1.0)*phi1 - (2.0*(b - c + d - e - 2.0)*phi1 - b + c - 2.0*d + 2.0*e + 4.0)*phi2 - b + c + 2.0)*phi3 + 2.0*c;
			double ratio = dnorm/norm2;
			x[0] = (-b*(phi2 - 1)*phi3 + ((phi2 - 1)*phi3 - phi2 + 1)*c + a*phi2)/norm - ratio*Q[0];
			x[1] = (-a*phi2)/norm - ratio*Q[1];
			x[2] = (b*(phi2 - 1)*phi3)/norm - ratio*Q[2];
			x[3] = (-((phi2 - 1)*phi3 - phi2 + 1)*c)/norm - ratio*Q[3];
			
			x[4] = a/norm - ratio*Q[4];
			x[5] = (-d*(phi2 - 1)*phi3 + ((phi2 - 1)*phi3 - phi2 + 1)*e - a)/norm - ratio*Q[5];
			x[6] = (d*(phi2 - 1)*phi3)/norm - ratio*Q[6];
			x[7] = (-((phi2 - 1)*phi3 - phi2 + 1)*e)/norm - ratio*Q[7];
			
			x[8] = b/norm - ratio*Q[8];
			x[9] = -d*phi2/norm - ratio*Q[9];
			x[10] = (d*phi2 + (phi2 - 1)*phi3 - b - phi2 + 1)/norm - ratio*Q[10];
			x[11] = (-(phi2 - 1)*phi3 + phi2 - 1)/norm - ratio*Q[11];
			
			x[12] = c/norm - ratio*Q[12];
			x[13] = -e*phi2/norm - ratio*Q[13];
			x[14] = (phi2 - 1)*phi3/norm - ratio*Q[14];
			x[15] = (e*phi2 - (phi2 - 1)*phi3 - c)/norm - ratio*Q[15];
		}
		// dQdphi2
		else if(index == 6){
			const double dnorm = -2*(a - c - e)*phi1*phi1 + 4*(phi1*phi1 - (phi1*phi1 - 2*phi1 + 1)*phi2 - 2*phi1 + 1)*phi3*phi3 + 2*(a - c - 2*e)*phi1 - 4*(e*phi1*phi1 - 2*e*phi1 + e)*phi2 + 2*((b - c + d - e - 2)*phi1*phi1 - (b - c + 2*d - 2*e - 4)*phi1 - 2*((d - e - 1)*phi1*phi1 - 2*(d - e - 1)*phi1 + d - e - 1)*phi2 + d - e - 2)*phi3 + 2*e;
			double ratio = dnorm/norm2;
			x[0] = (-b*(phi1 - 1)*phi3 + ((phi1 - 1)*phi3 - phi1 + 1)*c + a*(phi1 - 1))/norm - ratio*Q[0];
			x[1] = -a*(phi1 - 1)/norm - ratio*Q[1];
			x[2] = b*(phi1 - 1)*phi3/norm - ratio*Q[2];
			x[3] = -((phi1 - 1)*phi3 - phi1 + 1)*c/norm - ratio*Q[3];
			
			x[4] = -ratio*Q[4];
			x[5] = (-d*(phi1 - 1)*phi3 + ((phi1 - 1)*phi3 - phi1 + 1)*e)/norm - ratio*Q[5];
			x[6] = d*(phi1 - 1)*phi3/norm - ratio*Q[6];
			x[7] = -((phi1 - 1)*phi3 - phi1 + 1)*e/norm - ratio*Q[7];
			
			x[8] = -ratio*Q[8];
			x[9] = -d*(phi1 - 1)/norm - ratio*Q[9];
			x[10] = (d*(phi1 - 1) + (phi1 - 1)*phi3 - phi1 + 1)/norm - ratio*Q[10];
			x[11] = (-(phi1 - 1)*phi3 + phi1 - 1)/norm - ratio*Q[11];
			
			x[12] = -ratio*Q[12];
			x[13] = -e*(phi1 - 1)/norm - ratio*Q[13];
			x[14] = (phi1 - 1)*phi3/norm - ratio*Q[14];
			x[15] = (e*(phi1 - 1) - (phi1 - 1)*phi3)/norm - ratio*Q[15];
			
		}// dQdphi3
		else if(index == 7){
			const double dnorm = -2*(b - c - 1)*phi1*phi1 - 2*((d - e - 1)*phi1*phi1 - 2*(d - e - 1)*phi1 + d - e - 1)*phi2*phi2 + 2*(b - c - 2)*phi1 + 2*((b - c + d - e - 2)*phi1*phi1 - (b - c + 2*d - 2*e - 4)*phi1 + d - e - 2)*phi2 - 4*((phi1*phi1 - 2*phi1 + 1)*phi2*phi2 + phi1*phi1 - 2*(phi1*phi1 - 2*phi1 + 1)*phi2 - 2*phi1 + 1)*phi3 + 2;
			const double ratio = dnorm/norm2;
			x[0] = (-((phi1 - 1)*phi2 - phi1 + 1)*b + ((phi1 - 1)*phi2 - phi1 + 1)*c)/norm - ratio*Q[0];
			x[1] = - ratio*Q[1];
			x[2] = ((phi1 - 1)*phi2 - phi1 + 1)*b/norm - ratio*Q[2];
			x[3] = -((phi1 - 1)*phi2 - phi1 + 1)*c/norm - ratio*Q[3];
			
			x[4] = - ratio*Q[4];
			x[5] = (-((phi1 - 1)*phi2 - phi1 + 1)*d + ((phi1 - 1)*phi2 - phi1 + 1)*e)/norm - ratio*Q[5];
			x[6] = ((phi1 - 1)*phi2 - phi1 + 1)*d/norm - ratio*Q[6];
			x[7] = -((phi1 - 1)*phi2 - phi1 + 1)*e/norm - ratio*Q[7];
			
			x[8] = - ratio*Q[8];
			x[9] = - ratio*Q[9];
			x[10] = ((phi1 - 1)*phi2 - phi1 + 1)/norm - ratio*Q[10];
			x[11] = (-(phi1 - 1)*phi2 + phi1 - 1)/norm - ratio*Q[11];
			
			x[12] = - ratio*Q[12];
			x[13] = - ratio*Q[13];
			x[14] = ((phi1 - 1)*phi2 - phi1 + 1)/norm - ratio*Q[14];
			x[15] = (-(phi1 - 1)*phi2 + phi1 - 1)/norm - ratio*Q[15];
		}
		else{
			exit(1);
		}
		//m->dQ_need_update = false;
	}
	
	free(Q);
	
	const double *v = m->eigendcmp->eval;
	double xx[16];
	Matrix_mult3(mat, (const double**)m->eigendcmp->Invevec, x, 4,4,4,4);
	Matrix_mult4(xx, mat, (const double**)m->eigendcmp->evec, 4,4,4,4);
	// up to now the above operations can be recycled across branches
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			if(v[i] != v[j]){
				mat[i*4+j] = xx[i*4+j]*(exp(v[i]*t) - exp(v[j]*t))/(v[i]-v[j]);
			}
			else{
				mat[i*4+j] = xx[i*4+j]*t*exp(v[i]*t);
			}
		}
	}
	Matrix_mult3(xx, (const double**)m->eigendcmp->evec, mat, 4,4,4,4);
	Matrix_mult4(mat, xx, (const double**)m->eigendcmp->Invevec, 4,4,4,4);
}

