// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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
