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

#include "nucsubst.h"
#include "matrix.h"

static void _nuc_unrestricted_nonstat_update_Q( SubstitutionModel *m );

SubstitutionModel * new_NONSTATNucleotideModel(){
	Simplex* freqs = new_Simplex("nonstat", 4);
    SubstitutionModel *m = create_nucleotide_model("NONSTAT", NON_STATIONARY_DNA, freqs);
	
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
	
	
	m->update_Q = _nuc_unrestricted_nonstat_update_Q;
    return m;
}

SubstitutionModel * new_NONSTATNucleotideModel_with_parameters( Simplex* freqs, const Parameters *rates ){
	SubstitutionModel *m = create_nucleotide_model("NONSTAT", NON_STATIONARY_DNA, freqs);
	
	m->rates = new_Parameters( Parameters_count(rates) );
	for(int i = 0; i < Parameters_count(rates); i++){
		Parameters_add(m->rates, Parameters_at(rates, i) );
	}
	
	m->update_Q = _nuc_unrestricted_nonstat_update_Q;
	
	
	return m;
}

void _nuc_unrestricted_nonstat_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	const double* freqs = m->get_frequencies(m);
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
            m->Q[i][j] = m->Q[i][j] * freqs[j];
        }
    }
    
	make_zero_rows( m->Q, 4);
	normalize_Q( m->Q, freqs, 4 );
    m->need_update = false;
}
