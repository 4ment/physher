/*
 *  dayhoff.c
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

#include "dayhoff.h"

#include "matrix.h"

static void _dayoff_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	const double* f = m->get_frequencies(m);
	for ( int i = 0; i < m->nstate; i++ )  {
		for ( int j = i + 1; j < m->nstate; j++ ) {
			m->Q[i][j] = AMINO_ACID_MODEL_DAYHOFF[i][j] * f[j];
			m->Q[j][i] = AMINO_ACID_MODEL_DAYHOFF[i][j] * f[i];
		}
	}
	make_zero_rows( m->Q, 20);
	normalize_Q( m->Q, f, 20 );
	m->need_update = false;
}

SubstitutionModel * new_DAYHOFF(){
    return new_DAYHOFF_with_parameters(NULL);
}

SubstitutionModel *new_DAYHOFF_with_parameters(Parameter *freqs) {
    Parameter *freqs2 = NULL;
    if (freqs == NULL) {
        freqs2 = new_Parameter2("dayoff.freqs", AMINO_ACID_MODEL_DAYHOFF_FREQUENCIES,
                                20, new_Constraint(0.0, 1.0));
        Parameter_set_estimate(freqs2, false);
        freqs = freqs2;
    }
    SubstitutionModel *m = create_aa_model("DAYHOFF", DAYHOFF, freqs);
    if (freqs2 != NULL) {
        free_Parameter(freqs2);
    }
    m->update_Q = _dayoff_update_Q;
    _dayoff_update_Q(m);
    update_eigen_system(m);

    return m;
}
