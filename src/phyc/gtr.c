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

SubstitutionModel * new_GTR(){
    double freqs[4] = {0.25, 0.25, 0.25, 0.25};
    double rates[5] = {1,1,1,1,1};
    return new_GTR_with_values(freqs, rates);
}

SubstitutionModel * new_GTR_with_values( const double *freqs, const double *rates ){
    check_frequencies( freqs, 4 );
    
    Parameters* ratess = new_Parameters( 5 );
    Parameters_move(ratess, new_Parameter_with_postfix("gtr.a", "model", rates[0], new_Constraint(0.001, 100) ) );
    Parameters_move(ratess, new_Parameter_with_postfix("gtr.b", "model", rates[1], new_Constraint(0.001, 100) ) );
    Parameters_move(ratess, new_Parameter_with_postfix("gtr.c", "model", rates[2], new_Constraint(0.001, 100) ) );
    Parameters_move(ratess, new_Parameter_with_postfix("gtr.d", "model", rates[3], new_Constraint(0.001, 100) ) );
    Parameters_move(ratess, new_Parameter_with_postfix("gtr.e", "model", rates[4], new_Constraint(0.001, 100) ) );
    
    Parameters* freqss = new_Parameters( 3 );
    Parameters_move(freqss, new_Parameter_with_postfix("gtr.piA", "model", freqs[0]/freqs[3], new_Constraint(0.001, 100) ) );
    Parameters_move(freqss, new_Parameter_with_postfix("gtr.piC", "model", freqs[1]/freqs[3], new_Constraint(0.001, 100) ) );
    Parameters_move(freqss, new_Parameter_with_postfix("gtr.piG", "model", freqs[2]/freqs[3], new_Constraint(0.001, 100) ) );
	
	SubstitutionModel *m = new_GTR_with_parameters(freqss, ratess);
	free_Parameters(freqss);
	free_Parameters(ratess);
    return m;
}

SubstitutionModel * new_GTR_with_parameters( const Parameters* freqs, const Parameters* rates ){
	//check_frequencies( freqs, 4 );
	
	SubstitutionModel *m = create_nucleotide_model("GTR", GTR);
	
	m->_freqs = dvector(4);
	if(Parameters_count(freqs) == 4){
		m->update_frequencies = nucleotide_update_freqs;
		for (int i = 0; i < Parameters_count(freqs); i++) {
			m->_freqs[i] = Parameters_value(freqs, i);
		}
	}
	else{
		m->update_frequencies = nucleotide_update_freqs_relative;
		fprintf(stderr, "new_GTR_with_parameters need to be implemented\n");
		exit(1);
	}
	
	m->update_Q = gtr_update_Q;
	
	m->dPdp = gtr_dQ;
	m->dQ = dvector(16);
	
	m->rates = new_Parameters( 5 );
	for(int i = 0; i < Parameters_count(rates); i++){
		Parameters_add(m->rates, Parameters_at(rates, i) );
	}
	
	m->freqs = new_Parameters( 3 );
	for(int i = 0; i < Parameters_count(freqs); i++){
		Parameters_add(m->rates, Parameters_at(freqs, i) );
	}
	
	return m;
}


void gtr_update_Q( SubstitutionModel *m ){
	
    m->Q[0][1] = m->Q[1][0] = Parameters_value(m->rates, 0); // a
    m->Q[0][2] = m->Q[2][0] = Parameters_value(m->rates, 1); // b
    m->Q[0][3] = m->Q[3][0] = Parameters_value(m->rates, 2); // c
    
    m->Q[1][2] = m->Q[2][1] = Parameters_value(m->rates, 3); // d
    m->Q[1][3] = m->Q[3][1] = Parameters_value(m->rates, 4); // e
    
    m->Q[2][3] = m->Q[3][2] = 1; // f
    
    double temp;
    for ( int i = 0; i < m->nstate; i++ )  {
        for ( int j = i + 1; j < m->nstate; j++ ) {
            temp = m->Q[i][j];
            m->Q[i][j] = temp * m->_freqs[j];
            m->Q[j][i] = temp * m->_freqs[i];
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
        
        const double phi1 = Parameters_value(m->freqs, 0);
        const double phi2 = Parameters_value(m->freqs, 1);
        const double phi3 = Parameters_value(m->freqs, 2);

        const double piT = m->_freqs[3];
        
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

