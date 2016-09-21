/*
 *  nucsubst.c
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

#include "nucsubst.h"

#include "matrix.h"

static void nuc_sym_update_Q( SubstitutionModel *m );


SubstitutionModel * new_ReversibleNucleotideModel( const char model[5] ){
    SubstitutionModel *m = NULL;
    
    m = create_nucleotide_model("nucleotide", REVERSIBLE_DNA);
    m->_freqs = dvector(4);
    //default model has fixed frequencies
    for ( int i = 0; i < 4; i++ ) {
        m->_freqs[i] = 0.25;
    }
    
    int max = 0;
    for ( int i = 0; i < 5; i++ ) {
        int code = model[i];
        
        // lower case model e.g abaab
        if ( code >= 97 ) {
            code -= 97;
        }
        // uppper case model e.g ABAAB
        else if( code >= 65 ){
            code -= 65;
        }
        // numbers e.g 01001
        else if( code >= 48 ){
            code -= 48;
        }
        
        m->model[i] = code;
        
        if ( code > max ) {
            max = code;
        }
    }
    max++;
    
    if( max > 0 ){
        StringBuffer *buffer = new_StringBuffer(10);
        m->rates = new_Parameters( max );
        for ( int i = 0; i < max; i++ ) {
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "nuc.%d", i);
            Parameters_add(m->rates, new_Parameter_with_postfix(buffer->c, "model", 1, new_Constraint(0.001, 100) ) );
        }
        free_StringBuffer(buffer);
        
        m->update_Q = nuc_sym_update_Q;
    }
    
    return m;
}

void ReversibleNucleotideModel_estimate_freqs( SubstitutionModel *m ){
    if ( m->freqs == NULL ) {
        m->freqs = new_Parameters( 3 );
        double aux1 = m->_freqs[1] /   (1 - m->_freqs[0]);
        double aux2 = m->_freqs[2] / ( (1 - m->_freqs[0]) * (1 - aux1) );
        Parameters_add(m->freqs, new_Parameter_with_postfix("nuc.piA", "model", m->_freqs[0], new_Constraint(0.001, 0.999) ) );
        Parameters_add(m->freqs, new_Parameter_with_postfix("nuc.piC", "model", aux1,     new_Constraint(0.001, 0.999) ) );
        Parameters_add(m->freqs, new_Parameter_with_postfix("nuc.piT", "model", aux2,     new_Constraint(0.001, 0.999) ) );
        
        m->update_frequencies = nucleotide_update_freqs;
        m->need_update = true;
    }
}

void nuc_sym_update_Q( SubstitutionModel *m ){
    
    
    m->Q[0][1] = m->Q[1][0] = Parameters_value(m->rates, m->model[0]); // a
    m->Q[0][2] = m->Q[2][0] = Parameters_value(m->rates, m->model[1]); // b
    m->Q[0][3] = m->Q[3][0] = Parameters_value(m->rates, m->model[2]); // c
    
    m->Q[1][2] = m->Q[2][1] = Parameters_value(m->rates, m->model[3]); // d
    m->Q[1][3] = m->Q[3][1] = Parameters_value(m->rates, m->model[4]); // e
    
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

