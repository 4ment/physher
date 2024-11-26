/*
 *  mg94.c
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

#include "mg94.h"

#include <assert.h>

#include "geneticcode.h"
#include "matrix.h"

static void _mg_update_Q(SubstitutionModel* m);

SubstitutionModel* new_MG94(Parameter* freqs, unsigned gen_code) {
    return new_MG94_with_values(freqs, 1, 1, 1, gen_code);
}

SubstitutionModel* new_MG94_with_values(Parameter* freqs, const double alpha,
                                        const double beta, const double kappa,
                                        unsigned gen_code) {
    SubstitutionModel* m = create_codon_model("MG94", MG94, gen_code, freqs);

    // Functions
    m->update_Q = _mg_update_Q;

    m->rates = new_Parameters(3);
    Parameters_move(m->rates, new_Parameter_with_postfix("mg.kappa", "model", kappa,
                                                         new_Constraint(0.001, 100)));
    Parameters_move(m->rates, new_Parameter_with_postfix("mg.alpha", "model", alpha,
                                                         new_Constraint(0.001, 100)));
    Parameters_move(m->rates, new_Parameter_with_postfix("mg.beta", "model", beta,
                                                         new_Constraint(0.001, 100)));

    return m;
}

SubstitutionModel* new_MG94_with_parameters(Parameter* freqs, Parameter* alpha,
                                            Parameter* beta, Parameter* kappa,
                                            unsigned gen_code) {
    SubstitutionModel* m = create_codon_model("MG94", MG94, gen_code, freqs);

    // Functions
    m->update_Q = _mg_update_Q;

    m->rates = new_Parameters(3);
    Parameters_add(m->rates, kappa);
    Parameters_add(m->rates, alpha);
    Parameters_add(m->rates, beta);

    return m;
}

static void _mg_update_Q(SubstitutionModel* m) {
    if (!m->need_update) return;
    int i = 0;
    int j, jj;

    double alpha = Parameters_value(m->rates, 1);
    double beta = Parameters_value(m->rates, 2);
    double kappa_alpha = Parameters_value(m->rates, 0) * alpha;
    double kappa_beta = Parameters_value(m->rates, 0) * beta;
    const double* freqs = m->get_frequencies(m);

    for (int ii = 0; ii < 64; ii++) {
        if (GENETIC_CODE_TABLES[m->gen_code][ii] == '*') continue;
        m->Q[i][i] = 0.0;
        for (jj = ii + 1, j = i + 1; jj < 64; jj++) {
            if (GENETIC_CODE_TABLES[m->gen_code][jj] == '*') continue;

            // fprintf(stderr, "%d %d - %d %d\n", i,j,ii,jj);

            int subst = -1;
            m->Q[i][j] = m->Q[j][i] = 0.0;

            if (CODON_TRIPLETS[ii][0] != CODON_TRIPLETS[jj][0]) {
                subst = 0;
            }
            if (CODON_TRIPLETS[ii][1] != CODON_TRIPLETS[jj][1]) {
                if (subst != -1) {
                    j++;
                    continue;
                }
                subst = 1;
            }
            if (CODON_TRIPLETS[ii][2] != CODON_TRIPLETS[jj][2]) {
                if (subst != -1) {
                    j++;
                    continue;
                }
                subst = 2;
            }
            assert(subst != -1);  // there should be a substitution

            // there is one transition
            if ((CODON_TRIPLETS[ii][subst] == 'A' &&
                 CODON_TRIPLETS[jj][subst] == 'G')  // A G
                || (CODON_TRIPLETS[ii][subst] == 'G' &&
                    CODON_TRIPLETS[jj][subst] == 'A')  // G A
                || (CODON_TRIPLETS[ii][subst] == 'C' &&
                    CODON_TRIPLETS[jj][subst] == 'T')  // C T
                || (CODON_TRIPLETS[ii][subst] == 'T' &&
                    CODON_TRIPLETS[jj][subst] == 'C')) {  // T C

                // synonymous
                if (GENETIC_CODE_TABLES[m->gen_code][ii] ==
                    GENETIC_CODE_TABLES[m->gen_code][jj]) {
                    m->Q[i][j] = kappa_alpha * freqs[j];
                    m->Q[j][i] = kappa_alpha * freqs[i];
                } else {
                    m->Q[i][j] = kappa_beta * freqs[j];
                    m->Q[j][i] = kappa_beta * freqs[i];
                }
            }
            // there is one transversion
            else {
                // synonymous
                if (GENETIC_CODE_TABLES[m->gen_code][ii] ==
                    GENETIC_CODE_TABLES[m->gen_code][jj]) {
                    m->Q[i][j] = alpha * freqs[j];
                    m->Q[j][i] = alpha * freqs[i];
                } else {
                    m->Q[i][j] = beta * freqs[j];
                    m->Q[j][i] = beta * freqs[i];
                }
            }
            j++;
        }
        i++;
    }
    make_zero_rows(m->Q, m->nstate);
    normalize_Q(m->Q, freqs, m->nstate);
    m->need_update = false;
}
