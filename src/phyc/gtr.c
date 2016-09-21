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

SubstitutionModel * new_GTR_with_parameters( Parameters *freqs, Parameters *rates ){
    //check_frequencies( freqs, 4 );
    
    SubstitutionModel *m = create_nucleotide_model("GTR", GTR);
    
    m->_freqs = dvector(4);
    nucleotide_update_freqs(m);
    
    m->update_frequencies = nucleotide_update_freqs;
    
    m->rates = clone_Parameters(rates, true);
    m->freqs = clone_Parameters(freqs, true);
    return m;
}

SubstitutionModel * new_GTR(){
    double freqs[4] = {0.25, 0.25, 0.25, 0.25};
    double rates[5] = {1,1,1,1,1};
    return new_GTR_with_values(freqs, rates);
}

SubstitutionModel * new_GTR_with_values( const double *freqs, const double *rates ){
    check_frequencies( freqs, 4 );
    
    SubstitutionModel *m = create_nucleotide_model("GTR", GTR);
    
    m->_freqs = clone_dvector( freqs, 4);
    
    m->update_frequencies = nucleotide_update_freqs;
    m->update_Q = gtr_update_Q;
    
    m->rates = new_Parameters( 5 );
    Parameters_add(m->rates, new_Parameter_with_postfix("gtr.a", "model", rates[0], new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("gtr.b", "model", rates[1], new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("gtr.c", "model", rates[2], new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("gtr.d", "model", rates[3], new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("gtr.e", "model", rates[4], new_Constraint(0.001, 100) ) );
    
    m->freqs = new_Parameters( 3 );
    double aux1 = freqs[1] /   (1 - freqs[0]);
    double aux2 = freqs[2] / ( (1 - freqs[0]) * (1 - aux1) );
    Parameters_add(m->freqs, new_Parameter_with_postfix("gtr.piA", "model", freqs[0], new_Constraint(0.001, 0.999) ) );
    Parameters_add(m->freqs, new_Parameter_with_postfix("gtr.piC", "model", aux1,     new_Constraint(0.001, 0.999) ) );
    Parameters_add(m->freqs, new_Parameter_with_postfix("gtr.piT", "model", aux2,     new_Constraint(0.001, 0.999) ) );
    
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
//    Parameters_print(m->rates);
//    Parameters_print(m->freqs);
    update_eigen_system( m );
    m->need_update = false;
}

