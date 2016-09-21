/*
 *  nonstat.c
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

#include "nonstat.h"

#include "matrix.h"

static void _nuc_unrestricted_nonstat_update_Q( SubstitutionModel *m );

SubstitutionModel * new_NONSTATNucleotideModel(){
    double freqs[4] = {0.25,0.25,0.25,0.25};
    SubstitutionModel *m = new_NONSTATNucleotideModel_with_values(freqs);
    return m;
}

SubstitutionModel * new_NONSTATNucleotideModel_with_values( const double *freqs ){
    SubstitutionModel *m = NULL;
    
    m = create_nucleotide_model("NONSTAT", NON_STATIONARY_DNA);
    
    m->_freqs = clone_dvector(freqs, 4);
    
    m->freqs = new_Parameters( 3 );
    double aux1 = m->_freqs[1] /   (1 - m->_freqs[0]);
    double aux2 = m->_freqs[2] / ( (1 - m->_freqs[0]) * (1 - aux1) );
    Parameters_add(m->freqs, new_Parameter_with_postfix("gtr.piA", "model", m->_freqs[0], new_Constraint(0.001, 0.999) ) );
    Parameters_add(m->freqs, new_Parameter_with_postfix("gtr.piC", "model", aux1,     new_Constraint(0.001, 0.999) ) );
    Parameters_add(m->freqs, new_Parameter_with_postfix("gtr.piT", "model", aux2,     new_Constraint(0.001, 0.999) ) );
    
    m->rates = new_Parameters( 11 );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r1",  "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r2",  "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r3",  "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r4",  "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r5",  "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r6",  "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r7",  "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r8",  "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r9",  "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r10", "model", 1, new_Constraint(0.001, 100) ) );
    Parameters_add(m->rates, new_Parameter_with_postfix("unres.r11", "model", 1, new_Constraint(0.001, 100) ) );
    
    
    m->update_Q = _nuc_unrestricted_nonstat_update_Q;
    m->update_frequencies = nucleotide_update_freqs;
    
    
    return m;
}

void _nuc_unrestricted_nonstat_update_Q( SubstitutionModel *m ){
    
    int index = 0;
    for ( int i = 0; i < 4; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
            if( i != j ){
                if( index == 11 ) m->Q[i][j] = 1;
                else{
                    m->Q[i][j] = Parameters_value(m->rates, index);
                    index++;
                }
            }
        }
    }
    
    for ( int i = 0; i < m->nstate; i++ )  {
        for ( int j = 0; j < m->nstate; j++ ) {
            m->Q[i][j] = m->Q[i][j] * m->_freqs[j];
        }
    }
    
    update_eigen_system( m );
    m->need_update = false;
}
