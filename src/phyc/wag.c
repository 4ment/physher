/*
 *  wag.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 16/09/2014.
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

#include "wag.h"


#include "matrix.h"


SubstitutionModel * new_WAG(){
	Simplex* freqs = new_Simplex("wag", 20);
	freqs->set_values(freqs, AMINO_ACID_MODEL_WAG_FREQUENCIES);
    SubstitutionModel *m = new_WAG_with_parameters(freqs);
    
    return m;
}

SubstitutionModel * new_WAG_with_parameters( Simplex *freqs ){
	if(freqs == NULL){
		freqs = new_Simplex("wag", 20);
		freqs->set_values(freqs, AMINO_ACID_MODEL_WAG_FREQUENCIES);
	}
    SubstitutionModel *m = create_aa_model("WAG", WAG, freqs);
	
	const double* f = m->get_frequencies(m);
	
    for ( int i = 0; i < m->nstate; i++ )  {
        for ( int j = i + 1; j < m->nstate; j++ ) {
            m->Q[i][j] = AMINO_ACID_MODEL_WAG[i][j] * f[j];
            m->Q[j][i] = AMINO_ACID_MODEL_WAG[i][j] * f[i];
        }
    }
    
    update_eigen_system( m );
    m->need_update = false;
    
    return m;
}
