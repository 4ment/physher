// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "lg.h"

#include "matrix.h"

static void _lg_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	
	const double* f = m->get_frequencies(m);
	for ( int i = 0; i < m->nstate; i++ )  {
		for ( int j = i + 1; j < m->nstate; j++ ) {
			m->Q[i][j] = AMINO_ACID_MODEL_LG[i][j] * f[j];
			m->Q[j][i] = AMINO_ACID_MODEL_LG[i][j] * f[i];
		}
	}
	make_zero_rows( m->Q, 20);
	normalize_Q( m->Q, f, 20 );
	m->need_update = false;
}

SubstitutionModel * new_LG(){
    SubstitutionModel* m = new_LG_with_parameters(NULL);

    return m;
}

SubstitutionModel* new_LG_with_parameters(Parameter* freqs) {
    Parameter* freqs2 = NULL;
    if (freqs == NULL) {
        freqs2 = new_Parameter2("lg.freqs", AMINO_ACID_MODEL_LG_FREQUENCIES, 20,
                                new_Constraint(0.0, 1.0));
        Parameter_set_estimate(freqs, false);
        freqs = freqs2;
    }
    SubstitutionModel* m = create_aa_model("LG", LG, freqs);
    if (freqs2 != NULL) {
        free_Parameter(freqs2);
    }
    m->update_Q = _lg_update_Q;
    _lg_update_Q(m);
    update_eigen_system(m);

    return m;
}
