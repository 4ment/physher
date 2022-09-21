/*
 *  gtr.c
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

#include "gtr.h"

#include "nucsubst.h"
#include "matrix.h"

/***************
 Q matrix
 
 A C G T
 A  * a b c
 C  a * d e
 G  b d * 1
 T  c e 1 *
 
 freqs  A C G T
 0 1 2 3
 
 rates a b c d e
 0 1 2 3 4
  
 a = aa;
 b = (1 - aa) *      bb;
 c = (1 - aa) * (1 - bb) *      cc;
 d = (1 - aa) * (1 - bb) * (1 - cc);
 
 aa = a;
 bb = b /  (1 - aa);
 cc = c / ((1 - aa)*(1 - bb));
 
 ***************/

static void _gtr_dPdp(SubstitutionModel *m, int index, double* mat, double t);

/*! Initialize a general time reversible matrix.
 *  Return a Q matrix (pi+freq)
 */


SubstitutionModel * new_GTR(Simplex* freqs){
	Parameter* ac = new_Parameter_with_postfix("gtr.a", "model", 1, new_Constraint(0.001, 100) );
	Parameter* ag = new_Parameter_with_postfix("gtr.b", "model", 1, new_Constraint(0.001, 100) );
	Parameter* at = new_Parameter_with_postfix("gtr.c", "model", 1, new_Constraint(0.001, 100) );
	Parameter* cg = new_Parameter_with_postfix("gtr.d", "model", 1, new_Constraint(0.001, 100) );
	Parameter* ct = new_Parameter_with_postfix("gtr.e", "model", 1, new_Constraint(0.001, 100) );
	
	SubstitutionModel* model = new_GTR_with_parameters(freqs, ac, ag, at, cg, ct, NULL);
	
	free_Parameter(ac);
	free_Parameter(ag);
	free_Parameter(at);
	free_Parameter(cg);
	free_Parameter(ct);
	
	return model;
}

SubstitutionModel * new_GTR_with_values( const double *freqs, const double *rates ){
    check_frequencies( freqs, 4 );
    
	Parameter* ac = new_Parameter_with_postfix("gtr.a", "model", 1, new_Constraint(0.001, 100) );
	Parameter* ag = new_Parameter_with_postfix("gtr.b", "model", 1, new_Constraint(0.001, 100) );
	Parameter* at = new_Parameter_with_postfix("gtr.c", "model", 1, new_Constraint(0.001, 100) );
	Parameter* cg = new_Parameter_with_postfix("gtr.d", "model", 1, new_Constraint(0.001, 100) );
	Parameter* ct = new_Parameter_with_postfix("gtr.e", "model", 1, new_Constraint(0.001, 100) );
	
    Simplex* sfreqs = new_Simplex("gtr", 4);
	sfreqs->set_values(sfreqs, freqs);
	
	SubstitutionModel* model = new_GTR_with_parameters(sfreqs, ac, ag, at, cg, ct, NULL);

	free_Parameter(ac);
	free_Parameter(ag);
	free_Parameter(at);
	free_Parameter(cg);
	free_Parameter(ct);
    return model;
}

SubstitutionModel * new_GTR_with_parameters( Simplex* freqs, Parameter* ac, Parameter* ag, Parameter* at, Parameter* cg, Parameter* ct, Parameter* gt){
	
	SubstitutionModel *m = create_nucleotide_model("GTR", GTR, freqs);
	
	m->update_Q = gtr_update_Q;
	
	m->dPdp = _gtr_dPdp;
	m->dQ = dvector(16);
	m->relativeTo = -1;
	
	m->rates = new_Parameters(5);
	if (ac != NULL) {
		Parameters_add(m->rates, ac);
	}
	else{
		m->relativeTo = 0;
	}
	if (ag != NULL) {
		Parameters_add(m->rates, ag);
	}
	else{
		m->relativeTo = 1;
	}
	if (at != NULL) {
		Parameters_add(m->rates, at);
	}
	else{
		m->relativeTo = 2;
	}
	if (cg != NULL) {
		Parameters_add(m->rates, cg);
	}
	else{
		m->relativeTo = 3;
	}
	if (ct != NULL) {
		Parameters_add(m->rates, ct);
	}
	else{
		m->relativeTo = 4;
	}
	if (gt != NULL) {
		Parameters_add(m->rates, gt);
	}
	else{
		m->relativeTo = 5;
	}
	for(size_t i = 0; i < Parameters_count(m->rates); i++){
		Parameters_at(m->rates, i)->id = i;
	}
	return m;
}

SubstitutionModel * new_GTR_with_simplexes( Simplex* freqs, Simplex* rates){
	
	SubstitutionModel *m = create_nucleotide_model("GTR", GTR, freqs);
	
	m->update_Q = gtr_simplexes_update_Q;
	
	m->dPdp = _gtr_dPdp;
	m->dQ = dvector(16);
	
	m->rates_simplex = rates;
	
	return m;
}

void gtr_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	double temp;
	int index = 0;
	const double* freqs = m->get_frequencies(m);
	for ( int i = 0; i < 4; i++ )  {
		for ( int j = i + 1; j < 4; j++ ) {
			if(m->relativeTo != index){
				temp = Parameters_value(m->rates, index++);
			}
			else{
				temp = 1;
			}
			m->Q[i][j] = temp * freqs[j];
			m->Q[j][i] = temp * freqs[i];
		}
	}
	make_zero_rows( m->Q, 4);
	m->norm = normalize_Q( m->Q, freqs, 4 );
	m->need_update = false;
}

void gtr_simplexes_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	int index = 0;
	const double* freqs = m->simplex->get_values(m->simplex);
	const double* rates = m->rates_simplex->get_values(m->rates_simplex);
	for ( int i = 0; i < 4; i++ )  {
		for ( int j = i + 1; j < 4; j++ ) {
			m->Q[i][j] = rates[index] * freqs[j];
			m->Q[j][i] = rates[index] * freqs[i];
			index++;
		}
	}
	make_zero_rows( m->Q, 4);
	m->norm = normalize_Q( m->Q, freqs, 4 );
	m->need_update = false;
}

void nuc_sym_update_Q_unnormalized( SubstitutionModel *m, double* mat ){
	double temp;
	int index = 0;
	const double* freqs = m->get_frequencies(m);
	for ( int i = 0; i < 4; i++ )  {
		for ( int j = i + 1; j < 4; j++ ) {
			if(m->relativeTo != index){
				temp = Parameters_value(m->rates, index++);
			}
			else{
				temp = 1;
			}
			mat[i*4+j] = temp * freqs[j];
			mat[j*4+i] = temp * freqs[i];
		}
	}
	for (int i = 0; i < 4; i++) {
		double sum = 0.0;
		for (int j = 0; j < 4; j++) {
			if (i != j) {
				sum += mat[i*4+j];
			}
		}
		mat[i*4+i] = -sum;
	}
}

void gtr_setup_rates(double* Q, const double* rates, size_t dim){
	size_t index = 0;
	for (size_t i = 0; i < dim; i++) {
		Q[i*dim+i] = 0;
		for (size_t j = i + 1; j < dim; j++) {
			Q[i*dim+j] = Q[j*dim+i] = rates[index++];
		}
	}
}

void gtr_setup_rates_relative(double* Q, Parameters* rates, size_t dim, size_t relativeTo){
	size_t index = 0;
	double rate;
	for (size_t i = 0; i < dim; i++) {
		Q[i*dim+i] = 0;
		for (size_t j = i + 1; j < dim; j++) {
			if(relativeTo != index){
				rate = Parameters_value(rates, index++);
			}
			else{
				rate = 1;
			}
			Q[i*dim+j] = Q[j*dim+i] = rate;
		}
	}
}

static void _gtr_dQdp(SubstitutionModel *m, size_t index){
	double *dQ = m->dQ;
	size_t stateCount = m->nstate;
	memset(dQ, 0.0, sizeof(double)*stateCount*stateCount);
	if(m->need_update){
		m->update_Q(m);
	}

	const double* frequencies = m->simplex->get_values(m->simplex);
	double norm = m->norm;
	double dnorm = 0;
	if(m->rates_simplex != NULL && ((m->grad_wrt_reparam && index < m->rates_simplex->K - 1) || (!m->grad_wrt_reparam && index < m->rates_simplex->K))){
		double dR[6];
		if(m->grad_wrt_reparam){
			m->rates_simplex->gradient(m->rates_simplex, index, dR);
		}
		else{
			memset(dR, 0.0, sizeof(double)*6);
			dR[index] = 1.0;
		}
		gtr_setup_rates(dQ, dR, stateCount);
		build_Q_flat(dQ, frequencies, stateCount);
		dnorm = normalizing_constant_Q_flat(dQ, frequencies, stateCount);
	}
	else if(m->rates_simplex == NULL && index < Parameters_count(m->rates)){
		double dR[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		dR[index] = 1.0;
		gtr_setup_rates(dQ, dR, stateCount);
		build_Q_flat(dQ, frequencies, stateCount);
		dnorm = normalizing_constant_Q_flat(dQ, frequencies, stateCount);
	}
	else{
		size_t rateCount = -1;
		if(m->rates_simplex == NULL){
			rateCount = Parameters_count(m->rates);
		}
		else if(m->grad_wrt_reparam){
			rateCount = m->rates_simplex->K - 1;
		}
		else{
			rateCount = m->rates_simplex->K;
		}	
		if(m->rates_simplex != NULL){
			const double* rates = m->rates_simplex->get_values(m->rates_simplex);
			gtr_setup_rates(dQ, rates, stateCount);
		}
		else{
			gtr_setup_rates_relative(dQ, m->rates, stateCount, m->relativeTo);
		}
		double dF[4];
		if(m->grad_wrt_reparam){
			m->simplex->gradient(m->simplex, index - rateCount, dF);
		}
		else{
			memset(dF, 0.0, sizeof(double)*4);
			dF[index - rateCount] = 1.0;
		}
		build_Q_flat(dQ, dF, stateCount);
		
		for (size_t i = 0; i < stateCount; i++) {
			dnorm -= dQ[i*stateCount+i]*frequencies[i] + m->Q[i][i]*norm*dF[i];
		}
	}

	for(size_t i = 0; i < stateCount; i++){
		for(size_t j = 0; j < stateCount; j++){
			// (dQ . norm - norm_dQ . Q)/norm^2
			dQ[i*stateCount + j] = (dQ[i*stateCount + j] - m->Q[i][j]*dnorm)/norm;
		}
	}
}

static void _gtr_dPdp(SubstitutionModel *m, int index, double* mat, double t){
	if(m->need_update){
		m->update_Q(m);
		m->dQ_need_update = true;
	}
	if(m->dQ_need_update){
		_gtr_dQdp(m, index);
		m->dQ_need_update = false;
	}
	dPdp_with_dQdp(m, m->dQ, mat, t);
}
