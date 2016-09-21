/*
 *  lg.c
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

#include "lg.h"

#include "matrix.h"


SubstitutionModel * new_LG(){
    double freqs[20];
    for ( int i = 0; i < 20; i++ ) {
        freqs[i] = 0.05;
    }
    SubstitutionModel *m = new_LG_with_values(freqs);
    
    return m;
}

SubstitutionModel * new_LG_with_values( const double *freqs ){
    
    SubstitutionModel *m = create_aa_model("LG", LG);
    
    if( freqs != NULL ){
        check_frequencies( freqs, 20 );
        m->_freqs = clone_dvector( freqs, 20 );
    }
    else {
        m->_freqs = clone_dvector( AMINO_ACID_MODEL_LG_FREQUENCIES, 20 );
    }
    
    for ( int i = 0; i < m->nstate; i++ )  {
        for ( int j = i + 1; j < m->nstate; j++ ) {
            m->Q[i][j] = AMINO_ACID_MODEL_LG[i][j] * m->_freqs[j];
            m->Q[j][i] = AMINO_ACID_MODEL_LG[i][j] * m->_freqs[i];
        }
    }
    
    update_eigen_system( m );
    m->need_update = false;
    
    return m;
}
