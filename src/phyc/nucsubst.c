// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "nucsubst.h"

#include "matrix.h"

static void nuc_sym_update_Q( SubstitutionModel *m );

// we can give a string of 5 characters 00111 -> rates are relative to the last rate (gt)
// we can give a string of 6 characters 0*0111 -> rates are relative to the second rate
// use GTR if you want a simplex and 6 rates
SubstitutionModel * new_ReversibleNucleotideModel_with_parameters( const char* model, Parameter* freqs, const Parameters* rates){
    SubstitutionModel *m = NULL;
    
    m = create_nucleotide_model("nucleotide", REVERSIBLE_DNA, freqs);
	m->relativeTo = -1;
	m->rates = new_Parameters(5);
    
    int max = 0;
	m->model = new_DiscreteParameter(5);
	int index = 0;
    for ( int i = 0; i < strlen(model); i++ ) {
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
		// relative to this one
		else{
			m->relativeTo = i;
			continue;
		}
        
        m->model->values[index] = code;
        
        if ( code > max ) {
            max = code;
        }
		index++;
    }
    max++;
	
	if(strlen(model) == 5) m->relativeTo = 5;
    
    if( max > 0 ){
		Parameters_add_parameters(m->rates, rates);
		
        m->update_Q = nuc_sym_update_Q;
    }
    
    return m;
}

void nuc_sym_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	const unsigned* model = m->model->values;
	const double* freqs = Parameter_values(m->simplex);
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
    
    make_zero_rows( m->Q, 4);
	normalize_Q( m->Q, freqs, 4 );
    m->need_update = false;
}

