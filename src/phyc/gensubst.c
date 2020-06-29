/*
 *  gensubst.c
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

#include "gensubst.h"

#include "matrix.h"


// symmetric matrix
static void _p_t_integrated_exp( SubstitutionModel *m, const double t, double *P ){
	if( m->need_update ){
		double lambda = 10;
		double lambda2 = lambda*lambda;
		double p11 = 2.0/(9.*lambda2 + 36.*lambda + 32.); // 2 mutations
		double p01 = (3.*lambda + 2.)/(9.*lambda2 + 36.*lambda + 32.);// 1 mutation
		double p00 = (9.*lambda2 + 18.*lambda + 2.)/(9.*lambda2 + 36.*lambda + 32.);
		
		for (size_t i = 0; i < m->nstate; i++){
			m->PP[i][i] = p00;
			for (size_t j = i+1; j < m->nstate; j++){
				size_t subst = 0;
				if(m->datatype->states[i][0] != m->datatype->states[j][0]){
					subst++;
				}
				if(m->datatype->states[i][1] != m->datatype->states[j][1]){
					subst++;
				}
				if(subst == 1){
					m->PP[i][j] = m->PP[j][i] = p01;
				}
				else{
					m->PP[i][j] = m->PP[j][i] = p11;
				}
			}
		}
		m->need_update = false;
	}
	
    for ( size_t i = 0; i < m->nstate; i++){
		memcpy(P+i*m->nstate, m->PP[i], sizeof(double)*m->nstate);
	}
}

bool isSymmetric(const unsigned *matrix, size_t dim){
    for ( int i = 0; i < dim; i++ ) {
        for ( int j = i+1; j < dim; j++ ) {
            if( matrix[i*dim+j] != matrix[j*dim+i]){
                return false;
            }
        }
    }
    return true;
}


void _nonreversible_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	const double* freqs = m->simplex->get_values(m->simplex);
    int index = 0;
	const unsigned* model = m->model->values;
    for ( int i = 0; i < m->nstate; i++ )  {
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++]) * freqs[j];
        }
        j++;
        for ( ; j < m->nstate; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++]) * freqs[j];
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
	m->normalize = 1;
    if(m->normalize) m->norm = normalize_Q( m->Q, freqs, m->nstate );
    
    m->need_update = false;
}

void _nonreversible_update_Q2( SubstitutionModel *m ){
	if(!m->need_update) return;
	const double* freqs = m->simplex->get_values(m->simplex);
    int index = 0;
	const unsigned* model = m->model->values;
    for ( int i = 0; i < m->nstate; i++ )  {
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++]) * freqs[j];
        }
        j++;
        for ( ; j < m->nstate; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++]) * freqs[j];
        }
    }
    
    for ( int i = 0; i < m->nstate-1; i++ )  {
        m->Q[m->nstate-1][i] = 0;
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[m->nstate-1][i] += m->Q[i][j];
        }
        for ( j++ ; j < m->nstate; j++ ) {
            m->Q[m->nstate-1][i] += m->Q[i][j];
        }
        for ( j = 0; j < m->nstate-1; j++ ) {
            if( i != j ){
                m->Q[m->nstate-1][i] -= m->Q[j][i];
            }
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
	m->norm = 1;
    if(m->normalize)m->norm = normalize_Q( m->Q, freqs, m->nstate );
    
    m->need_update = false;
}


void _general_jc69_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	
	double offdiag = 1.0/(m->nstate-1);
	for ( int i = 0; i < m->nstate; i++ )  {
		m->Q[i][i] = -1;
		for ( int j = i + 1; j < m->nstate; j++ ) {
			m->Q[i][j] = m->Q[j][i] = offdiag;
		}
	}
    
    m->need_update = false;
}

void _general_jc69_p_t( SubstitutionModel *m, const double t, double *P ){
	double nstate = m->nstate;
    double temp = exp(-nstate/(nstate-1.0) * t);
	double p00 = 1.0/nstate + (nstate-1.0)/nstate * temp;
	double p01 = 1.0/nstate - temp/nstate;
	for ( int i = 0; i < m->nstate; i++ )  {
		P[i*m->nstate+i] = p00;
		for ( int j = i + 1; j < m->nstate; j++ ) {
			P[i*m->nstate+j] = P[j*m->nstate+i] = p01;
		}
	}
}

void _general_jc69_dp_dt( SubstitutionModel *m, const double t, double *P ){
    double nstate = m->nstate;
    double temp = exp(-nstate/(nstate-1.0) * t);
	for ( int i = 0; i < m->nstate; i++ )  {
		P[i*m->nstate+i] = -temp;
		for ( int j = i + 1; j < m->nstate; j++ ) {
			P[i*m->nstate+j] = P[j*m->nstate+i] = temp/(nstate-1.0);
		}
	}
}

void _general_jc69_d2p_dt2( SubstitutionModel *m, const double t, double *P ){
    double nstate = m->nstate;
    double temp = exp(-nstate/(nstate-1.0) * t);
	for ( int i = 0; i < m->nstate; i++ )  {
		P[i*m->nstate+i] = temp * nstate/(nstate-1);
		for ( int j = i + 1; j < m->nstate; j++ ) {
			P[i*m->nstate+j] = P[j*m->nstate+i] = -temp*nstate/((nstate-1)*(nstate-1));
		}
	}
}

void _reversible_update_Q( SubstitutionModel *m ){
	if(!m->need_update) return;
	const unsigned* model = m->model->values;
	double temp;
	int index = 0;
	m->normalize = 1;
	
	if(m->simplex != NULL){
		const double* freqs = m->simplex->get_values(m->simplex);
		for ( int i = 0; i < m->nstate; i++ )  {
			index += i+1;
			for ( int j = i + 1; j < m->nstate; j++ ) {
				if(m->relativeTo != model[index]){
					temp = Parameters_value(m->rates, model[index]);
				}
				else{
					temp = 1;
				}
				index++;
				m->Q[i][j] = temp * freqs[j];
				m->Q[j][i] = temp * freqs[i];
			}
		}
		make_zero_rows( m->Q, m->nstate);
		if(m->normalize)
			m->normalize = normalize_Q( m->Q, freqs, m->nstate );
	}
	else{
		for ( int i = 0; i < m->nstate; i++ )  {
			index += i+1;
			for ( int j = i + 1; j < m->nstate; j++ ) {
				if(m->relativeTo != model[index]){
					temp = Parameters_value(m->rates, model[index]);
				}
				else{
					temp = 1;
				}
				index++;
				m->Q[i][j] = temp;
				m->Q[j][i] = temp;
			}
		}
		make_zero_rows( m->Q, m->nstate);
	}
    
    m->need_update = false;
}

SubstitutionModel * new_GeneralModel_with_parameters( DataType* datatype, DiscreteParameter* model, const Parameters* rates, Simplex* freqs, int relativeTo, bool normalize ){
    //bool sym = isSymmetric(model->values, freqs->K);
		size_t dim = sqrt(model->length);
	bool sym = isSymmetric(model->values, dim);
    SubstitutionModel *m = NULL;
    if(sym){
        m = create_general_model("GENERAL", REVERSIBLE, datatype, freqs);
        m->update_Q = _reversible_update_Q;
    }
	else{
        m = create_general_model("UREVGENERAL", NONREVERSIBLE, datatype, freqs);
        m->update_Q = _nonreversible_update_Q;
    }
    m->normalize = normalize;
	m->relativeTo = relativeTo;
    
	m->model = model;
	
	m->rates = new_Parameters(Parameters_count(rates));
	Parameters_add_parameters(m->rates, rates);
    
    return m;
}


SubstitutionModel * new_GeneralJC69Model_with_parameters( DataType* datatype, Simplex* freqs_simplex, bool normalize ){
    SubstitutionModel *m = create_general_model("GENERAL", REVERSIBLE, datatype, freqs_simplex);
	m->update_Q = _general_jc69_update_Q;
	m->p_t = m->p_t_transpose = _general_jc69_p_t;
	m->dp_dt = m->dp_dt_transpose = _general_jc69_dp_dt;
	m->d2p_d2t = m->d2p_d2t_transpose = _general_jc69_d2p_dt2;
//	m->p_t = m->p_t_transpose = _p_t_integrated_exp;
    m->normalize = normalize;
	m->update_Q(m);
    return m;
}
