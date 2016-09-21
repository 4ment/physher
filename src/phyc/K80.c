/*
 *  K80.c
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

#include "K80.h"

#include "matrix.h"

/***************
 K80
 Q matrix
 
    A C G T
 A  * 1 1 1
 C  1 * 1 k
 G  k 1 * 1
 T  1 k 1 *
 
 freqs  A C G T = 0.25
 0 1 2 3
 
 rates a b c d e
 0 1 2 3 4
 ***************/

static void k80_update_Q( SubstitutionModel *m );

SubstitutionModel * new_K80(){
    SubstitutionModel *m = create_nucleotide_model("K80", K80);
    
    m->_freqs = dvector(4);
    for ( int i = 0; i < 4; i++ ) m->_freqs[i] = 0.25;
    
    m->update_Q = k80_update_Q;
    
    m->rates = new_Parameters( 1 );
    Parameters_add(m->rates, new_Parameter_with_postfix("k80.kappa", "model", 3, new_Constraint(0.0001, 100) ) );
    
    return m;
}

SubstitutionModel * new_K80_with_values( const double kappa ){
    SubstitutionModel *m = create_nucleotide_model("K80", K80);
    
    m->_freqs = dvector(4);
    for ( int i = 0; i < 4; i++ ) m->_freqs[i] = 0.25;
    
    m->update_Q = k80_update_Q;
    
    m->rates = new_Parameters( 1 );
    Parameters_add(m->rates, new_Parameter_with_postfix("k80.kappa", "model", kappa, new_Constraint(0.0001, 100) ) );
    
    return m;
}

void k80_update_Q( SubstitutionModel *m ){
    
    m->Q[1][3] = m->Q[3][1] = m->Q[0][2] = m->Q[2][0] = Parameters_value(m->rates, 0)*0.25; // kappa
    m->Q[0][1] = m->Q[1][0] = m->Q[0][3] = m->Q[3][0] = m->Q[1][2] = m->Q[2][1] = m->Q[2][3] = m->Q[3][2] = 0.25;
    
    make_zero_rows( m->Q, 4);
    //normalize_Q( m->Q, m->_freqs, m->nstate );
    double subst = -(m->Q[0][0]*0.25 + m->Q[1][1]*0.25 + m->Q[2][2]*0.25 + m->Q[3][3]*0.25);
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            m->Q[i][j] /= subst;
        }
    }
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    m->need_update = false;
}
