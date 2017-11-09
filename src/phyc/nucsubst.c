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


SubstitutionModel * new_ReversibleNucleotideModel( const char model[5], Simplex* freqs ){
    SubstitutionModel *m = NULL;
    
    m = create_nucleotide_model("nucleotide", REVERSIBLE_DNA, freqs);
	m->relativeTo = 5;
    
    int max = 0;
	m->model = new_DiscreteParameter("nucleotide.model", 5);
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
        
        m->model->values[i] = code;
        
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
            Parameters_move(m->rates, new_Parameter_with_postfix(buffer->c, "model", 1, new_Constraint(0.001, 100) ) );
        }
        free_StringBuffer(buffer);
        
        m->update_Q = nuc_sym_update_Q;
    }
    
    return m;
}

void nuc_sym_update_Q( SubstitutionModel *m ){
	const unsigned* model = m->model->values;
	const double* freqs = m->simplex->get_values(m->simplex);
    double temp;
	int index = 0;
    for ( int i = 0; i < 4; i++ )  {
        for ( int j = i + 1; j < 4; j++ ) {
			if(m->relativeTo != index){
				temp = Parameters_value(m->rates, model[index++]);
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

