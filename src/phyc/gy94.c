/*
 *  gy94.c
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

#include "gy94.h"

#include <assert.h>

#include "matrix.h"
#include "geneticcode.h"

static void _gy_update_Q( SubstitutionModel *m );

SubstitutionModel *new_GY94(Parameter *freqs, unsigned gen_code) {
    SubstitutionModel *m = new_GY94_with_values(freqs, 1, 1, gen_code);
    return m;
}

SubstitutionModel *new_GY94_with_values(Parameter *freqs, const double omega,
                                        const double kappa, unsigned gen_code) {
    SubstitutionModel *m = create_codon_model("GY94", GY94, gen_code, freqs);

    // Functions
    m->update_Q = _gy_update_Q;

    m->rates = new_Parameters(2);
    Parameters_move(m->rates, new_Parameter_with_postfix("gy.kappa", "model", kappa,
                                                         new_Constraint(0.001, 100)));
    Parameters_move(m->rates, new_Parameter_with_postfix("gy.omega", "model", omega,
                                                         new_Constraint(0.001, 100)));

    return m;
}

static void _gy_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
    int i = 0;
    int j,jj;
    
    double kappa = Parameters_value(m->rates, 0);
    double omega = Parameters_value(m->rates, 1);
	const double* freqs = m->get_frequencies(m);
    
    for ( int ii = 0; ii < 64; ii++ ) {
        if ( GENETIC_CODE_TABLES[m->gen_code][ii] == '*' ) continue;
        m->Q[i][i] = 0.0;
        for ( jj = ii+1, j=i+1; jj < 64; jj++ ) {
            if ( GENETIC_CODE_TABLES[m->gen_code][jj] == '*' ) continue;
            
            //fprintf(stderr, "%d %d - %d %d\n", i,j,ii,jj);
            
            int subst = -1;
            m->Q[i][j] = m->Q[j][i] = 0.0;
            
            if ( CODON_TRIPLETS[ii][0] != CODON_TRIPLETS[jj][0] ) {
                subst = 0;
            }
            if ( CODON_TRIPLETS[ii][1] != CODON_TRIPLETS[jj][1] ) {
                if( subst != -1 ){
                    j++;
                    continue;
                }
                subst = 1;
            }
            if ( CODON_TRIPLETS[ii][2] != CODON_TRIPLETS[jj][2] ) {
                if( subst != -1 ){
                    j++;
                    continue;
                }
                subst = 2;
            }
            assert( subst != -1 ); // there should be a substitution
            
            // there is one transition
            if (   (CODON_TRIPLETS[ii][subst] == 'A' && CODON_TRIPLETS[jj][subst] == 'G')     // A G
                || (CODON_TRIPLETS[ii][subst] == 'G' && CODON_TRIPLETS[jj][subst] == 'A')     // G A
                || (CODON_TRIPLETS[ii][subst] == 'C' && CODON_TRIPLETS[jj][subst] == 'T')     // C T
                || (CODON_TRIPLETS[ii][subst] == 'T' && CODON_TRIPLETS[jj][subst] == 'C') ) { // T C
                
                // synonymous
                if ( GENETIC_CODE_TABLES[m->gen_code][ii] == GENETIC_CODE_TABLES[m->gen_code][jj] ) {
                    m->Q[i][j] = kappa * freqs[j];
                    m->Q[j][i] = kappa * freqs[i];
                }
                else {
                    m->Q[i][j] = kappa * omega * freqs[j];
                    m->Q[j][i] = kappa * omega * freqs[i];
                }
            }
            // there is one transversion
            else{
                // synonymous
                if ( GENETIC_CODE_TABLES[m->gen_code][ii] == GENETIC_CODE_TABLES[m->gen_code][jj] ) {
                    m->Q[i][j] = freqs[j];
                    m->Q[j][i] = freqs[i];
                }
                else {
                    m->Q[i][j] = omega * freqs[j];
                    m->Q[j][i] = omega * freqs[i];
                }
            }
            j++;
        }
        i++;
    }
	
    make_zero_rows( m->Q, m->nstate);
	normalize_Q( m->Q, freqs, m->nstate );
    
    m->need_update = false;
}



void print_GY( unsigned gen_code ){
    int i = 0;
    int j,jj;
    
    int **Q = imatrix(NUMBER_OF_CODONS[gen_code], NUMBER_OF_CODONS[gen_code]);
    
    for ( int ii = 0; ii < 64; ii++ ) {
        if ( GENETIC_CODE_TABLES[gen_code][ii] == '*' ) continue;
        Q[i][i] = -1;
        
        for ( jj = ii+1, j=i+1; jj < 64; jj++ ) {
            if ( GENETIC_CODE_TABLES[gen_code][jj] == '*' ) continue;
            
            //fprintf(stderr, "%d %d - %d %d\n", i,j,ii,jj);
            
            int subst = -1;
            Q[i][j] = Q[j][i] = 0;
            
            if ( CODON_TRIPLETS[ii][0] != CODON_TRIPLETS[jj][0] ) {
                subst = 0;
            }
            if ( CODON_TRIPLETS[ii][1] != CODON_TRIPLETS[jj][1] ) {
                if( subst != -1 ){
                    j++;
                    continue;
                }
                subst = 1;
            }
            if ( CODON_TRIPLETS[ii][2] != CODON_TRIPLETS[jj][2] ) {
                if( subst != -1 ){
                    j++;
                    continue;
                }
                subst = 2;
            }
            assert( subst != -1 ); // there should be a substitution
            
            // there is one transition
            if ( (CODON_TRIPLETS[ii][subst] == 'A' && CODON_TRIPLETS[jj][subst] == 'G')       // A G
                || (CODON_TRIPLETS[ii][subst] == 'G' && CODON_TRIPLETS[jj][subst] == 'A')     // G A
                || (CODON_TRIPLETS[ii][subst] == 'C' && CODON_TRIPLETS[jj][subst] == 'T')     // C T
                || (CODON_TRIPLETS[ii][subst] == 'T' && CODON_TRIPLETS[jj][subst] == 'C') ) { // T C
                
                // synonymous
                if ( GENETIC_CODE_TABLES[gen_code][ii] == GENETIC_CODE_TABLES[gen_code][jj] ) {
                    Q[i][j] = 1;
                    Q[j][i] = 1;
                }
                else {
                    Q[i][j] = 3;
                    Q[j][i] = 3;
                }
            }
            // there is one transversion
            else{
                // synonymous
                if ( GENETIC_CODE_TABLES[gen_code][ii] == GENETIC_CODE_TABLES[gen_code][jj] ) {
                    Q[i][j] = 2;
                    Q[j][i] = 2;
                }
                else {
                    Q[i][j] = 4;
                    Q[j][i] = 4;
                }
            }
            assert(j < NUMBER_OF_CODONS[gen_code]);
            j++;
        }
        assert(i < NUMBER_OF_CODONS[gen_code]);
        i++;
    }
    //fprintf(stderr, "gy_update_Q %d %d\n",i,j);
    
    for ( int i = 0; i < 64; i++ ) {
        if ( GENETIC_CODE_TABLES[gen_code][i] == '*' ){
            continue;
        }
        fprintf(stdout, "\t%s", CODON_TRIPLETS[i]);
    }
    fprintf(stdout, "\n");
    for ( int i = 0; i < 64; i++ ) {
        if ( GENETIC_CODE_TABLES[gen_code][i] == '*' ){
            continue;
        }
        fprintf(stdout, "\t%c  ", GENETIC_CODE_TABLES[gen_code][i]);
    }
    
    fprintf(stdout, "\n");
    
    for ( int i = 0; i < NUMBER_OF_CODONS[gen_code]; i++ ) {
        fprintf(stdout, "%s\t", CODON_TRIPLETS[i]);
        for ( int j = 0; j < NUMBER_OF_CODONS[gen_code]; j++ ) {
            fprintf(stdout, "%d\t", Q[i][j]);
        }
        fprintf(stdout, "\n");
    }
    fflush(stdout);
    for ( int i = 0; i < NUMBER_OF_CODONS[gen_code]; i++ ) {
        free(Q[i]);
    }
    free(Q);
    printf("print_gy94\n");
    exit(2);
}
