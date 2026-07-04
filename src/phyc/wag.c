// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "wag.h"


#include "matrix.h"

static void _wag_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	
	const double* f = m->get_frequencies(m);
	for ( int i = 0; i < m->nstate; i++ )  {
		for ( int j = i + 1; j < m->nstate; j++ ) {
			m->Q[i][j] = AMINO_ACID_MODEL_WAG[i][j] * f[j];
			m->Q[j][i] = AMINO_ACID_MODEL_WAG[i][j] * f[i];
		}
	}
	make_zero_rows( m->Q, 20);
	normalize_Q( m->Q, f, 20 );
	m->need_update = false;
}

SubstitutionModel *new_WAG() { return new_WAG_with_parameters(NULL); }

SubstitutionModel *new_WAG_with_parameters(Parameter *freqs) {
    Parameter *freqs2 = NULL;
    if (freqs == NULL) {
        freqs2 = new_Parameter2("wag.freqs", AMINO_ACID_MODEL_WAG_FREQUENCIES, 20,
                                new_Constraint(0.0, 1.0));
        Parameter_set_estimate(freqs2, false);
    }
    SubstitutionModel *m = create_aa_model("WAG", WAG, freqs);
    m->update_Q = _wag_update_Q;
    _wag_update_Q(m);
    update_eigen_system(m);
    if (freqs2 != NULL) {
        free_Parameter(freqs2);
    }

    return m;
}
