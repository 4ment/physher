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
	
    Simplex* sfreqs = new_Simplex(4);
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


void gtr_dQ(SubstitutionModel *m, int index, double* mat, double t){
	
    double *x = m->dQ;
    
    if(m->need_update){
        m->update_Q(m);
    }
    
    if(m->dQ_need_update){
        const double a = Parameters_value(m->rates, 0);
        const double b = Parameters_value(m->rates, 1);
        const double c = Parameters_value(m->rates, 2);
        const double d = Parameters_value(m->rates, 3);
        const double e = Parameters_value(m->rates, 4);
        
        const double phi1 = Parameters_value(m->simplex->parameters, 0);
        const double phi2 = Parameters_value(m->simplex->parameters, 1);
        const double phi3 = Parameters_value(m->simplex->parameters, 2);

        const double piT = m->get_frequencies(m)[3];
        
        const double norm = piT*(phi1*(a*phi2 + b*phi3 + c) + phi2*(a*phi1 + d*phi3 + e) + phi3*(b*phi1 + d*phi2 + 1.0) + (c*phi1 + e*phi2 + phi3));
        const double norm2 = norm*norm;
        
        memset(x, 0, sizeof(double)*16);
        // derivative of norm with respect to a
        if(index == 0){
            const double dnorm = piT*2.0*phi1*phi2;
            const double ratio = dnorm/norm2;
            x[0] = (-norm*phi2 + (a*phi2 + b*phi3 + c)*dnorm)/norm2;
            x[1] = (norm*phi2 - a*phi2*dnorm)/norm2;
            x[2] = -b*phi3*ratio;
            x[3] = -c*ratio;
            
            x[4] = (norm*phi1 - a*phi1*dnorm)/norm2;
            x[5] = (-norm*phi1 + (a*phi1 + d*phi3 + e)*dnorm)/norm2;
            x[6] = -d*phi3*ratio;
            x[7] = -e*ratio;
            
            x[8] = -b*phi1*ratio;
            x[9] = -d*phi2*ratio;
            x[10] = (b*phi1 + d*phi2 + 1.0)*ratio;
            x[11] = -ratio;
            
            x[12] = -c*phi1*ratio;
            x[13] = -e*phi2*ratio;
            x[14] = -phi3*ratio;
            x[15] = (c*phi1 + e*phi2 + phi3)*ratio;
        }
        // derivative of norm with respect to b
        else if(index == 1){
            const double dnorm = piT*2.0*phi1*phi3;
            const double ratio = dnorm/norm2;
            x[0] = (-norm*phi3 + (a*phi2 + b*phi3 + c)*dnorm)/norm2;
            x[1] = -a*phi2*ratio;
            x[2] = (norm*phi3 - b*phi3*dnorm)/norm2;
            x[3] = -c*ratio;
            
            x[4] = -a*phi1*ratio;
            x[5] = ((a*phi1 + d*phi3 + e)*dnorm)/norm2;
            x[6] = -d*phi3*ratio;
            x[7] = -e*ratio;
            
            x[8] = (norm*phi1 - b*phi1*dnorm)/norm2;
            x[9] = -d*phi2*ratio;
            x[10] = (-norm*phi1 + (b*phi1 + d*phi2 + 1.0)*dnorm)/norm2;
            x[11] = -ratio;
            
            x[12] = -c*phi1*ratio;
            x[13] = -e*phi2*ratio;
            x[14] = -phi3*ratio;
            x[15] = (c*phi1 + e*phi2 + phi3)*ratio;
        }
        // derivative of norm with respect to c
        else if(index == 2){
            const double dnorm = piT*2.0*phi1;
            const double ratio = dnorm/norm2;
            x[0] = (-norm + (a*phi2 + b*phi3 + c)*dnorm)/norm2;
            x[1] = -a*phi2*ratio;
            x[2] = -b*phi3*ratio;
            x[3] = (norm - c*dnorm)/norm2;
            
            x[4] = -a*phi1*ratio;
            x[5] = (a*phi1 + d*phi3 + e)*dnorm/norm2;
            x[6] = -d*phi3*ratio;
            x[7] = -e*ratio;
            
            x[8] = -b*phi1*ratio;
            x[9] = -d*phi2*ratio;
            x[10] = (b*phi1 + d*phi2 + 1.0)*ratio;
            x[11] = -ratio;
            
            x[12] = (norm*phi1 - (c*phi1)*dnorm)/norm2;
            x[13] = -e*phi2*ratio;
            x[14] = -phi3*ratio;
            x[15] = (-norm*phi1 + (c*phi1 + e*phi2 + phi3)*dnorm)/norm2;
        }
        // derivative of norm with respect to d
        else if(index == 3){
            const double dnorm = piT*2.0*phi3*phi2;
            const double ratio = dnorm/norm2;
            x[0] = (a*phi2 + b*phi3 + c)*ratio;
            x[1] = -a*phi2*ratio;
            x[2] = -b*phi3*ratio;
            x[3] = -c*ratio;
            
            x[4] = -a*phi1*ratio;
            x[5] = (-norm*phi3 + (a*phi1 + d*phi3 + e)*dnorm)/norm2;
            x[6] = (norm*phi3 -d*phi3*dnorm)/norm2;
            x[7] = -e*ratio;
            
            x[8] = -b*phi1*ratio;
            x[9] = (norm*phi2 -d*phi2*dnorm)/norm2;
            x[10] = (-norm*phi2 +(b*phi1 + d*phi2 + 1.0)*dnorm)/norm2;
            x[11] = -ratio;
            
            x[12] = -c*phi1*ratio;
            x[13] = -e*phi2*ratio;
            x[14] = -phi3*ratio;
            x[15] = (c*phi1 + e*phi2 + phi3)*ratio;
        }
        // derivative of norm with respect to e
        else if(index == 4){
            const double dnorm = piT*2.0*phi2;
            const double ratio = dnorm/norm2;
            x[0] = (a*phi2 + b*phi3 + c)*ratio;
            x[1] = -a*phi2*ratio;
            x[2] = -b*phi3*ratio;
            x[3] = -c*ratio;
            
            x[4] = -a*phi1*ratio;
            x[5] = (-norm + (a*phi1 + d*phi3 + e)*dnorm)/norm2;
            x[6] = -d*phi3*ratio;
            x[7] = (norm - e*dnorm)/norm2;
            
            x[8] = -b*phi1*ratio;
            x[9] = -d*phi2*ratio;
            x[10] = (b*phi1 + d*phi2 + 1.0)*ratio;
            x[11] = -ratio;
            
            x[12] = -c*phi1*ratio;
            x[13] = (norm*phi2 - e*phi2*dnorm)/norm2;
            x[14] = -phi3*ratio;
            x[15] = (-norm*phi2 + (c*phi1 + e*phi2 + phi3)*dnorm)/norm2;
        }
        // dQdphi1
        else if(index == 5){
            const double dnorm = piT*(2.0*(a*phi2 + b*phi3 + c) - norm);
            const double ratio = dnorm/norm2;
            x[0] = (a*phi2 + b*phi3 + c)*ratio;
            x[1] = -a*phi2*ratio;
            x[2] = -b*phi3*ratio;
            x[3] = -c*ratio;
            
            x[4] = (norm*a - a*phi1*dnorm)/norm2;
            x[5] = (-norm*a + (a*phi1 + d*phi3 + e)*dnorm)/norm2;
            x[6] = -d*phi3*ratio;
            x[7] = -e*ratio;
            
            x[8] = (norm*b - b*phi1*dnorm)/norm2;
            x[9] = -d*phi2*ratio;
            x[10] = (-norm*b + (b*phi1 + d*phi2 + 1)*dnorm)/norm2;
            x[11] = -ratio;
            
            x[12] = (norm*c - c*phi1*dnorm)/norm2;
            x[13] = -e*phi2*ratio;
            x[14] = -phi3*ratio;
            x[15] = (-norm*c + (c*phi1 + e*phi2 + phi3)*dnorm)/norm2;
        }// dQdphi2
        else if(index == 6){
            const double dnorm = piT*(2.0*(a*phi1 + d*phi3 + e) - norm);
            const double ratio = dnorm/norm2;
            x[0] = (-norm*a + (a*phi2 + b*phi3 + c)*dnorm)/norm2;
            x[1] = (norm*a - a*phi2*dnorm)/norm2;
            x[2] = -b*phi3*ratio;
            x[3] = -c*ratio;
            
            x[4] = -a*phi1*ratio;
            x[5] = (a*phi1 + d*phi3 + e)*ratio;
            x[6] = -d*phi3*ratio;
            x[7] = -e*ratio;
            
            x[8] = -b*phi1*ratio;
            x[9] = (norm*d - d*phi2*dnorm)/norm2;
            x[10] = (-norm*d + (b*phi1 + d*phi2 + 1.0)*dnorm)/norm2;
            x[11] = -ratio;
            
            x[12] = -c*phi1*ratio;
            x[13] = (norm*e - e*phi2*dnorm)/norm2;
            x[14] = -phi3*ratio;
            x[15] = (-norm*e + (c*phi1 + e*phi2 + phi3)*dnorm)/norm2;
        }// dQdphi3
        else if(index == 7){
            const double dnorm = piT*(2.0*(b*phi1 + d*phi2 + 1.0) - norm);
            const double ratio = dnorm/norm2;
            x[0] = (-norm*b + (a*phi2 + b*phi3 + c)*dnorm)/norm2;
            x[1] = -a*phi2*ratio;
            x[2] = (norm*b - b*phi3*dnorm)/norm2;
            x[3] = -c*ratio;
            
            x[4] = -a*phi1*ratio;
            x[5] = (-norm*d + (a*phi1 + d*phi3 + e)*dnorm)/norm2;
            x[6] = (norm*d - d*phi3*dnorm)/norm2;
            x[7] = -e*ratio;
            
            x[8] = -b*phi1*ratio;
            x[9] = -d*phi2*ratio;
            x[10] = (b*phi1 + d*phi2 + 1)*ratio;
            x[11] = -ratio;
            
            x[12] = -c*phi1*ratio;
            x[13] = -e*phi2*ratio;
            x[14] = (norm - phi3*dnorm)/norm2;
            x[15] = (-norm + (c*phi1 + e*phi2 + phi3)*dnorm)/norm2;
        }
        else{
            exit(1);
        }
        //m->dQ_need_update = false;
    }
    
    
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

