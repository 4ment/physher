/*
 *  jc69.c
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

#include "jc69.h"

#include "matrix.h"

static void jc69_p_t( SubstitutionModel *m, const double t, double *P );
static double jc69_pij_t( SubstitutionModel *m, const int i, const int j, const double t );

static void jc69_dp_dt( SubstitutionModel *m, const double t, double *P );
static void jc69_d2p_dt2( SubstitutionModel *m, const double t, double *P );

SubstitutionModel * new_JC69(Simplex* freqs){
	for (int i = 0; i < Parameters_count(freqs->parameters); i++) {
		Parameters_set_estimate(freqs->parameters, false, i);
	}
    SubstitutionModel *m = create_substitution_model("JC69", JC69, DATA_TYPE_NUCLEOTIDE, freqs);
    m->nstate = 4;
    m->pij_t = jc69_pij_t;
    m->p_t   = m->p_t_transpose = jc69_p_t;
    m->dp_dt = m->dp_dt_transpose = jc69_dp_dt;
    m->d2p_d2t = m->d2p_d2t_transpose = jc69_d2p_dt2;
    return m;
}

double jc69_pij_t( SubstitutionModel *m, const int i, const int j, const double t ){
    double temp = exp(- 4./3. * t)*0.25;
    //	if(i == j) return 1./4. + 3./4. * exp(- 4./3. * t);
    //	else       return 1./4. - 1./4. * exp(- 4./3. * t);
    if(i == j) return 1./4. + 3*temp;
    else       return 1./4. - temp;
}

void jc69_p_t_rolled( SubstitutionModel *m, const double t, double *P ){
    int i,j;
    int k = 0;
    for ( i = 0; i < 4; i++) {
        for ( j = 0; j < 4; j++){
            if( i == j ) P[k] = 1./4. + 3./4. * exp(- 4./3. * t);
            else P[k] = 1./4. - 1./4. * exp(- 4./3. * t);
            k++;
        }
    }
}

void jc69_p_t( SubstitutionModel *m, const double t, double *P ){
    double temp = exp(- 4./3. * t);
    //    P[0] = P[5] = P[10] = P[15] = 1./4. + 3./4. * exp(- 4./3. * t);
    //    P[1] = P[2] = P[3] = P[4] = P[6] = P[7] = P[8] = P[9] = P[11] = P[12] = P[13] = P[14] = 1./4. - 1./4. * exp(- 4./3. * t);
    P[0] = P[5] = P[10] = P[15] = 0.25 + 3./4 * temp;
    P[1] = P[2] = P[3] = P[4] = P[6] = P[7] = P[8] = P[9] = P[11] = P[12] = P[13] = P[14] = 0.25 - 0.25 * temp;
}

void jc69_dp_dt( SubstitutionModel *m, const double t, double *P ){
    double temp = exp(- 4./3. * t);
    //    P[0] = P[5] = P[10] = P[15] = 1./4. + 3./4. * exp(- 4./3. * t);
    //    P[1] = P[2] = P[3] = P[4] = P[6] = P[7] = P[8] = P[9] = P[11] = P[12] = P[13] = P[14] = 1./4. - 1./4. * exp(- 4./3. * t);
    P[0] = P[5] = P[10] = P[15] = -temp;
    P[1] = P[2] = P[3] = P[4] = P[6] = P[7] = P[8] = P[9] = P[11] = P[12] = P[13] = P[14] = temp/3.0;
}

void jc69_d2p_dt2( SubstitutionModel *m, const double t, double *P ){
    double temp = exp(- 4./3. * t);
    P[0] = P[5] = P[10] = P[15] = temp * 4.0/3.0;
    P[1] = P[2] = P[3] = P[4] = P[6] = P[7] = P[8] = P[9] = P[11] = P[12] = P[13] = P[14] = -temp*4.0/9.0;
}

StringBuffer * SubstitutionModel_JC69_bufferize( StringBuffer *buffer, SubstitutionModel *m ){
	StringBuffer_append_string(buffer, "{");
	StringBuffer_append_strings(buffer, 2, "name:\"", m->name, "\"\n");
	StringBuffer_append_string(buffer, "}");
	
	return buffer;
}
