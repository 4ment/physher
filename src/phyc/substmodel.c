/*
 *  substmodel.m
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/4/10.
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


#include "substmodel.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <float.h>

#include "eigen.h"
#include "utils.h"
#include "matrix.h"
#include "mstring.h"
#include "geneticcode.h"

#include "model.h"

#include "wag.h"
#include "dayhoff.h"
#include "lg.h"

#include "hashtable.h"

#pragma mark Public function definition

#ifdef LISTENERS

static void _substitution_model_handle_change( Model *self, Model *model, int index ){
	SubstitutionModel* m = (SubstitutionModel*)self->obj;
	m->need_update = true;
	//TODO: should be smarter
	m->update_frequencies(m);
	m->dQ_need_update = true;
	self->listeners->fire( self->listeners, self, index );
}

static void _substitution_model_free( Model *self ){
	if(self->ref_count == 1){
		printf("Free subsitution model %s\n", self->name);
		SubstitutionModel* m = (SubstitutionModel*)self->obj;
	//	free_SubstitutionModel(m);
		free(m->name);
		if( m->_freqs != NULL ) free(m->_freqs);
		if( m->eigendcmp != NULL ) free_EigenDecomposition(m->eigendcmp);
		if( m->Q != NULL ) free_dmatrix(m->Q, m->nstate);
		if( m->PP != NULL ) free_dmatrix(m->PP, m->nstate);
		free_Parameters(m->rates);
		free_Parameters(m->freqs);
		if( m->model != NULL ) free(m->model);
		if( m->dQ != NULL ) free(m->dQ);
		free(m);
		
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _substitution_model_clone(Model* self, Hashtable *hash){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	SubstitutionModel* subst = (SubstitutionModel*)self->obj;
	SubstitutionModel* mclone = clone_substitution_model(subst);
	return new_SubstitutionModel2(self->name, mclone);
}

// SubstitutionModel2 listen to the rate and freq parameters
Model * new_SubstitutionModel2( const char* name, SubstitutionModel *sm ){
	Model *model = new_Model(name, sm);
	int i = 0;
	if ( sm->rates != NULL ) {
		for ( i = 0; i < Parameters_count(sm->rates); i++ ) {
			Parameters_at(sm->rates, i)->listeners->add( Parameters_at(sm->rates, i)->listeners, model );
		}
		Parameters_add_parameters(model->parameters, sm->rates);
	}
	if ( sm->freqs != NULL ) {
		for ( i = 0; i < Parameters_count(sm->freqs); i++ ) {
			Parameters_at(sm->freqs, i)->listeners->add( Parameters_at(sm->freqs, i)->listeners, model );
		}
		Parameters_add_parameters(model->parameters, sm->freqs);
	}
	
	model->update = _substitution_model_handle_change;
	model->free = _substitution_model_free;
	model->clone = _substitution_model_clone;
	return model;
}
#endif


double get_frequency( SubstitutionModel *m, int base ){
	return m->_freqs[base];
}

double * get_frequencies( SubstitutionModel *m ){
	return m->_freqs;
}


#pragma mark -
#pragma mark Private functions

static void foo_update( SubstitutionModel *m ){} // does not nothing

void general_update_freqs( SubstitutionModel *model ){
    model->need_update = true;
    model->dQ_need_update = true;
    model->_freqs[model->nstate-1] = 0;
    for(int i = 0; i < Parameters_count(model->freqs); i++){
        model->_freqs[model->nstate-1] += Parameters_value(model->freqs, i);
    }
    model->_freqs[model->nstate-1] = 1.0/(1.0+model->_freqs[model->nstate-1]);
    for(int i = 0; i < Parameters_count(model->freqs); i++){
        model->_freqs[i] = Parameters_value(model->freqs, i) * model->_freqs[model->nstate-1];
    }
}

static void _set_frequencies( SubstitutionModel *m, const double *freqs ){
	memcpy(m->_freqs, freqs, sizeof(double)*m->nstate);
	Parameters_set_value( m->freqs, 0, freqs[0]/freqs[m->nstate-1] ); // fire changes only once
	for (int i = 1; i < m->nstate-1; i++) {
		Parameters_at(m->freqs, i)->value = freqs[i]/freqs[m->nstate-1];
	}
	m->need_update = true;
	m->dQ_need_update = true;
}

static void _set_rates( SubstitutionModel *m, const double *rates ){
	for ( int i = 0; i < Parameters_count(m->rates); i++) {
		Parameters_set_value(m->rates,i, rates[i]);
	}
	m->need_update = true;
    m->dQ_need_update = true;
}

// Set the relative frequency (ie. Parameter)!!
static void _set_relative_frequency( SubstitutionModel *m, const double freq, const int index ){
	Parameters_set_value(m->freqs,index, freq);
	m->update_frequencies(m);
	m->need_update = true;
    m->dQ_need_update = true;
}

static void _set_rate( SubstitutionModel *m, const double rate, const int index ){
	Parameters_set_value(m->rates,index, rate);
	m->need_update = true;
    m->dQ_need_update = true;
}

// row major
static void _p_t( SubstitutionModel *m, const double t, double *P ){
    
    if( m->need_update ){
        m->update_Q(m);
    }
    
    int i,j,k;
    double *pP = P;
    
    double temp;
    
    if( isnan(t) ){
        for ( i = 0; i < m->nstate; i++){
            temp = P[m->nstate*m->nstate+i];//exp(m->eigendcmp->eval[i] * t);
            for ( j = 0; j < m->nstate; j++){
                m->PP[i][j] = m->eigendcmp->Invevec[i][j] * temp;
            }
        }
    }
    else {
        for ( i = 0; i < m->nstate; i++){
            temp = exp(m->eigendcmp->eval[i] * t);
            for ( j = 0; j < m->nstate; j++){
                m->PP[i][j] = m->eigendcmp->Invevec[i][j] * temp;
            }
        }
    }
    
    for ( i = 0; i < m->nstate; i++){
        for ( j = 0; j < m->nstate; j++){
            *pP = 0.;
            for ( k = 0; k < m->nstate; k++ )
                *pP += m->PP[k][j] * m->eigendcmp->evec[i][k];
            *pP = fabs(*pP);
            pP++;
        }
    }
    
}

// t(P D P-1) = t(P-1) D t(P)
// or column major
static void _p_t_transpose( SubstitutionModel *m, const double t, double *P ){
    
    if( m->need_update ){
        m->update_Q(m);
    }
    
    int i,j,k,l;
    l = 0;
    
    double temp;
    
    if( isnan(t) ){
        for ( i = 0; i < m->nstate; i++){
            temp = P[m->nstate*m->nstate+i];//exp(m->eigendcmp->eval[i] * t);
            for ( j = 0; j < m->nstate; j++){
                m->PP[i][j] = m->eigendcmp->evec[j][i] * temp;
            }
        }
    }
    else {
        for ( i = 0; i < m->nstate; i++){
            temp = exp(m->eigendcmp->eval[i] * t);
            for ( j = 0; j < m->nstate; j++){
                m->PP[i][j] = m->eigendcmp->evec[j][i] * temp;
            }
        }
    }
    
    for ( i = 0; i < m->nstate; i++){
        for ( j = 0; j < m->nstate; j++){
            P[l] = 0.;
            for ( k = 0; k < m->nstate; k++ ){
                P[l] += m->eigendcmp->Invevec[k][i] * m->PP[k][j];
            } 
            P[l] = fabs(P[l]);
            l++;
        }
    }
    
    
}

// row major
static void _aligned_p_t( SubstitutionModel *m, const double t, double *P ){
    
    if( m->need_update ){
        m->update_Q(m);
    }
    
    int i,j,k;
    double *pP = P;
    
    double temp;
    
    if( isnan(t) ){
        for ( i = 0; i < m->nstate; i++){
            temp = P[m->nstate*m->nstate+i];//exp(m->eigendcmp->eval[i] * t);
            for ( j = 0; j < m->nstate; j++){
                m->PP[i][j] = m->eigendcmp->Invevec[i][j] * temp;
            }
        }
    }
    else {
        for ( i = 0; i < m->nstate; i++){
            temp = exp(m->eigendcmp->eval[i] * t);
            for ( j = 0; j < m->nstate; j++){
                m->PP[i][j] = m->eigendcmp->Invevec[i][j] * temp;
            }
        }
    }
    
    for ( i = 0; i < m->nstate; i++){
        for ( j = 0; j < m->nstate; j++){
            *pP = 0.;
            for ( k = 0; k < m->nstate; k++ )
                *pP += m->PP[k][j] * m->eigendcmp->evec[i][k];
            *pP = fabs(*pP);
            pP++;
        }
        pP++;
    }
    
}

// t(P D P-1) = t(P-1) D t(P)
// or column major
static void _aligned_p_t_transpose( SubstitutionModel *m, const double t, double *P ){
    
    if( m->need_update ){
        m->update_Q(m);
    }
    
    int i,j,k,l;
    l = 0;
    
    double temp;
    
    if( isnan(t) ){
        for ( i = 0; i < m->nstate; i++){
            temp = P[m->nstate*m->nstate+i];//exp(m->eigendcmp->eval[i] * t);
            for ( j = 0; j < m->nstate; j++){
                m->PP[i][j] = m->eigendcmp->evec[j][i] * temp;
            }
        }
    }
    else {
        for ( i = 0; i < m->nstate; i++){
            temp = exp(m->eigendcmp->eval[i] * t);
            for ( j = 0; j < m->nstate; j++){
                m->PP[i][j] = m->eigendcmp->evec[j][i] * temp;
            }
        }
    }
    
    for ( i = 0; i < m->nstate; i++){
        for ( j = 0; j < m->nstate; j++){
            P[l] = 0.;
            for ( k = 0; k < m->nstate; k++ ){
                P[l] += m->eigendcmp->Invevec[k][i] * m->PP[k][j];
            } 
            P[l] = fabs(P[l]);
            l++;
        }
        l++;
    }
    
    
}

// P'(t)=QP(t) = DΛD-1 De^ΛtD-1 = DΛe^ΛtD-1
// row major
static void _dp_dt( SubstitutionModel *m, const double t, double *P ){
	
	if( m->need_update ){
		m->update_Q(m);
	}
	
	int i,j,k;
    double *pP = P;
	
	double temp;
	for ( i = 0; i < m->nstate; i++){
		temp = m->eigendcmp->eval[i] * exp(m->eigendcmp->eval[i] * t);
		for ( j = 0; j < m->nstate; j++){
			m->PP[i][j] = m->eigendcmp->Invevec[i][j] * temp;
		}
	}
	
	for ( i = 0; i < m->nstate; i++){
		for ( j = 0; j < m->nstate; j++){
			*pP = 0.;
			for ( k = 0; k < m->nstate; k++ ){
				*pP += m->PP[k][j] * m->eigendcmp->evec[i][k];
            }
            pP++;
		}
	}
	
}

// P'(t)=QP(t) = Q De^ΛtD-1
static void _dp_dt_v2( SubstitutionModel *m, const double t, double *P ){
	
	if( m->need_update ){
		m->update_Q(m);
	}
	
	int i,j,k;
    double *pP = P;
	
	double temp;
	for ( i = 0; i < m->nstate; i++){
		temp = exp(m->eigendcmp->eval[i] * t);
		for ( j = 0; j < m->nstate; j++){
			m->PP[i][j] = m->eigendcmp->Invevec[i][j] * temp;
		}
	}
    
	double *temp2 = dvector(m->nstate*m->nstate);
	for ( i = 0; i < m->nstate; i++){
		for ( j = 0; j < m->nstate; j++){
            *temp2 = 0;
			for ( k = 0; k < m->nstate; k++ ){
                *temp2 += m->PP[k][j] * m->eigendcmp->evec[i][k];
            }
            temp2++;
		}
	}
    temp2 -= 16;
	
	for( i = 0; i < m->nstate; i++ ){
		for( j = 0; j < m->nstate; j++ ){
            *pP = 0.;
			for( k = 0; k < m->nstate; k++ ){
                *pP += m->Q[k][j] * temp2[i*m->nstate+k];
            }
            pP++;
		}
	}
    free(temp2);
}

// t(P D P-1) = t(P-1) D t(P)
// or column major
static void _dp_dt_transpose( SubstitutionModel *m, const double t, double *P ){
	
	if( m->need_update ){
		m->update_Q(m);
	}
	
	int i,j,k,l;
	l = 0;

	double temp;
	for ( i = 0; i < m->nstate; i++){
		temp = m->eigendcmp->eval[i] * exp(m->eigendcmp->eval[i] * t);
		for ( j = 0; j < m->nstate; j++){
			m->PP[i][j] = m->eigendcmp->evec[j][i] * temp;
		}
	}
	
	for ( i = 0; i < m->nstate; i++){
		for ( j = 0; j < m->nstate; j++){
			P[l] = 0.;
			for ( k = 0; k < m->nstate; k++ ){
				P[l] += m->eigendcmp->Invevec[k][i] * m->PP[k][j];
			}
			l++;
		}
	}
	
}

// row major
static void _d2p_d2t( SubstitutionModel *m, const double t, double *P ){
	
	if( m->need_update ){
		m->update_Q(m);
	}
	
	int i,j,k;
    double *pP = P;
	
	double temp;
	for ( i = 0; i < m->nstate; i++){
		temp = m->eigendcmp->eval[i] * m->eigendcmp->eval[i] * exp(m->eigendcmp->eval[i] * t);
		for ( j = 0; j < m->nstate; j++){
			m->PP[i][j] = m->eigendcmp->Invevec[i][j] * temp;
		}
	}
	
	for ( i = 0; i < m->nstate; i++){
		for ( j = 0; j < m->nstate; j++){
			*pP = 0.;
			for ( k = 0; k < m->nstate; k++ )
				*pP += m->PP[k][j] * m->eigendcmp->evec[i][k];
			pP++;
		}
	}
	
}

// t(P D P-1) = t(P-1) D t(P)
// or column major
static void _d2p_d2t_transpose( SubstitutionModel *m, const double t, double *P ){
	
	if( m->need_update ){
		m->update_Q(m);
	}
	
	int i,j,k,l;
	l = 0;
	
	double temp;
	for ( i = 0; i < m->nstate; i++){
		temp = m->eigendcmp->eval[i] * m->eigendcmp->eval[i] * exp(m->eigendcmp->eval[i] * t);
		for ( j = 0; j < m->nstate; j++){
			m->PP[i][j] = m->eigendcmp->evec[j][i] * temp;
		}
	}
	
	for ( i = 0; i < m->nstate; i++){
		for ( j = 0; j < m->nstate; j++){
			P[l] = 0.;
			for ( k = 0; k < m->nstate; k++ ){
				P[l] += m->eigendcmp->Invevec[k][i] * m->PP[k][j];
			}
			l++;
		}
	}
	
}

static double _pij_t( SubstitutionModel *m, const int i, const int j, const double t ){
	int k;
	double pij = 0;
	
	if( m->need_update ){
		m->update_Q (m);
	}
	
	for ( k = 0; k < m->nstate; k++)
		pij += exp( m->eigendcmp->eval[k] * t ) * m->eigendcmp->evec[i][k] * m->eigendcmp->Invevec[k][j];
	
	if( pij < DBL_MIN ) return DBL_MIN;
	return pij;
}

#pragma mark -


SubstitutionModel * create_substitution_model( const char *name, const modeltype modelname, const datatype dtype ){
    SubstitutionModel *m = (SubstitutionModel *)malloc( sizeof(SubstitutionModel) );
    assert(m);
    m->id = 0;
    
    m->nstate = 0;
    m->name = String_clone(name);
    m->modeltype = modelname;
    m->dtype = dtype;
    
    // Parameters and Q matrix
    m->_freqs = NULL;
    m->Q = NULL;
    m->PP = NULL;
    
    m->rates = NULL;
    m->freqs = NULL;
    
    m->eigendcmp = NULL;
    
    m->need_update = true;
    m->dQ_need_update = true;
    
	// Functions
	m->set_frequencies = _set_frequencies; // real frequencies that sum to 1
	m->set_rates = _set_rates;
	
	m->set_relative_frequency = _set_relative_frequency; // relative frequency with dimention nstate-1
    m->set_rate = _set_rate;
    m->update_frequencies = foo_update;
    m->update_Q = foo_update;
    
    m->pij_t = _pij_t;
    m->p_t   = _p_t;
    m->p_t_transpose   = _p_t_transpose;
    m->dp_dt   = _dp_dt;
    m->dp_dt_transpose   = _dp_dt_transpose;
    
    m->d2p_d2t   = _d2p_d2t;
    m->d2p_d2t_transpose   = _d2p_d2t_transpose;
    
    m->dPdp  = NULL;
    m->dQ = NULL;
    
    m->free = free_SubstitutionModel;
    
    m->clone = clone_substitution_model;
    
    m->gen_code = 0;
    m->reversible = true;
    m->normalize = true;
    
    m->model = NULL;
    
    return m;
}

SubstitutionModel * create_nucleotide_model( const char *name, const modeltype modelname ){
    SubstitutionModel *m = create_substitution_model(name, modelname, DATA_TYPE_NUCLEOTIDE);
    m->nstate = 4;
    m->Q = dmatrix(m->nstate, m->nstate);
    m->PP = dmatrix(m->nstate, m->nstate);
    m->eigendcmp = new_EigenDecomposition(m->nstate);
    m->model = uivector(5);
    return m;
}

SubstitutionModel * create_codon_model( const char *name, const modeltype modelname, unsigned gen_code ){
    SubstitutionModel *m = create_substitution_model(name, modelname, DATA_TYPE_CODON);
    m->nstate = NUMBER_OF_CODONS[gen_code];
    m->gen_code = gen_code;
    m->Q = dmatrix(m->nstate, m->nstate);
    m->PP = dmatrix(m->nstate, m->nstate);
    m->eigendcmp = new_EigenDecomposition(m->nstate);
    return m;
}

SubstitutionModel * create_aa_model( const char *name, const modeltype modelname ){
    SubstitutionModel *m = create_substitution_model(name, modelname, DATA_TYPE_AMINO_ACID);
    m->nstate = 20;
    m->Q = dmatrix(m->nstate, m->nstate);
    m->PP = dmatrix(m->nstate, m->nstate);
    m->eigendcmp = new_EigenDecomposition(m->nstate);
    return m;
}

SubstitutionModel * create_general_model( const char *name, const modeltype modelname, int dim ){
    SubstitutionModel *m = create_substitution_model(name, modelname, DATA_TYPE_GENERIC);
    m->nstate = dim;
    m->Q = dmatrix(m->nstate, m->nstate);
    m->PP = dmatrix(m->nstate, m->nstate);
    m->eigendcmp = new_EigenDecomposition(m->nstate);
    return m;
}

void free_SubstitutionModel( SubstitutionModel *m){
	free_SubstitutionModel_share(m, false, false);
}

void free_SubstitutionModel_share( SubstitutionModel *m, bool share_freqs, bool share_rates ){
	free(m->name);
	if( m->_freqs != NULL ) free(m->_freqs);
	if( m->eigendcmp != NULL ) free_EigenDecomposition(m->eigendcmp);
	if( m->Q != NULL ) free_dmatrix(m->Q, m->nstate);
	if( m->PP != NULL ) free_dmatrix(m->PP, m->nstate);
	if( m->rates != NULL && !share_rates ) free_Parameters(m->rates);
	if( m->freqs != NULL && !share_freqs ) free_Parameters(m->freqs);
    if( m->model != NULL ) free(m->model);
    if( m->dQ != NULL ) free(m->dQ);
	free(m);
}


SubstitutionModel * clone_substitution_model(SubstitutionModel *m){
    return clone_substitution_model_share(m, false, false);
}

SubstitutionModel * clone_substitution_model_share(SubstitutionModel *m, bool share_freqs, bool share_rates){
    SubstitutionModel *clone =  create_substitution_model( m->name, m->modeltype, m->dtype);
    clone->nstate = m->nstate;
    clone->_freqs = clone_dvector(m->_freqs, m->nstate);
    
    if( m->Q != NULL ){
        clone->eigendcmp = clone_EigenDecomposition(m->eigendcmp);
        clone->Q = clone_dmatrix( m->Q, m->nstate, m->nstate );
        clone->PP = clone_dmatrix( m->PP, m->nstate, m->nstate );
    }
    
    clone->gen_code = m->gen_code;
    clone->reversible = m->reversible;
    clone->normalize = m->normalize;
    
    // Functions
    clone->pij_t = m->pij_t;
    clone->p_t   = m->p_t;
    clone->p_t_transpose   = m->p_t_transpose;
    
    clone->dp_dt   = m->dp_dt;
    clone->dp_dt_transpose   = m->dp_dt_transpose;
    
    clone->d2p_d2t   = m->d2p_d2t;
    clone->d2p_d2t_transpose   = m->d2p_d2t_transpose;
    
    m->dPdp = m->dPdp;
    if( m->dQ != NULL ) clone->dQ = clone_dvector(clone->dQ, m->nstate*m->nstate);
    
    clone->update_frequencies = m->update_frequencies;
	clone->update_Q = m->update_Q;
	clone->set_frequencies = m->set_frequencies;
	clone->set_rates = m->set_rates;
	clone->set_relative_frequency = m->set_relative_frequency;
	clone->set_rate = m->set_rate;
	
    clone->free = m->free;
    clone->clone = m->clone;
    
    if( m->rates != NULL ) clone->rates = clone_Parameters(m->rates, true );
    if( m->freqs != NULL ) clone->freqs = clone_Parameters(m->freqs, true );
    
    if( m->rates != NULL ){
        if( share_rates ){
            clone->rates = m->rates;
        }
        else {
            clone->rates = clone_Parameters(m->rates, true );
        }
    }
    if( m->freqs != NULL ){
        if( share_freqs ){
            clone->freqs = m->freqs;
        }
        else {
            clone->freqs = clone_Parameters(m->freqs, true );
        }
    }
    
    if( m->model != NULL ) clone->model = clone_uivector(m->model, m->nstate*m->nstate);
    
    clone->need_update = false;
    clone->dQ_need_update = false;
    
    return clone;
}


/*char  * SubstitutionModel_stringify( SubstitutionModel *m ){
	StringBuffer *buffer = new_StringBuffer(500);

	buffer = SubstitutionModel_bufferize( buffer, m );
	
	char *final = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	fprintf(stderr, "%d\n", (int)strlen(final));
	return final;
}*/

StringBuffer * SubstitutionModel_bufferize( StringBuffer *buffer, SubstitutionModel *m ){
	int i = 0;
	
	StringBuffer_append_string(buffer, "{");
	StringBuffer_append_strings(buffer, 2, "id:\"", m->id, "\",\n");
	StringBuffer_append_strings(buffer, 2, "name:\"", m->name, "\",\n");
	
	if(Parameters_count(m->rates) > 0){
		StringBuffer_append_string(buffer, "rates:[");
		for ( i = 0; i < Parameters_count(m->rates); i++) {
			StringBuffer_append_format(buffer," %f", Parameters_value(m->rates, i) );
		}
		StringBuffer_append_string(buffer, "]\n");
	}
	if(Parameters_count(m->rates) > 0){
		StringBuffer_append_string(buffer, "freqs:[");
		for ( i = 0; i < m->nstate; i++) {
			buffer = StringBuffer_append_format(buffer," %f", m->_freqs[i]);		
		}
		StringBuffer_append_string(buffer, "]\n");
	}
	StringBuffer_append_string(buffer, "}");

	return buffer;
}

/*void * SubstitutionModel_SML_to_object( ObjectStore *store, SMLNode node ){
	fprintf(stderr, "Model_SML_to_object\n");
	SubstitutionModel *m = NULL;
	
	char *rates_string = SML_get_data_of_child( node, "rates");
	char *freqs_string = SML_get_data_of_child( node, "freqs");
	
	char *type     = SML_get_data_of_child( node, "type");
	
	if( strcmp( type, "GTR") == 0 ){
		double *freqs = sml_string_to_double_array( freqs_string, 4 );
		double *rates = sml_string_to_double_array( rates_string, 5 );
		
		m = new_GTR_with_values( freqs, rates);
	}
	else if( strcmp( type, "HKY") == 0 ){
		double *kappa = sml_string_to_double_array( rates_string, 1 );
		double *freqs = sml_string_to_double_array( freqs_string, 4 );
		
		m = new_HKY_with_values(freqs, *kappa);
	}
	else if( strcmp( type, "JC69") ==0 ){
			m = new_JC69();
	}
	else if( strcmp( type, "CUSTOM") ==0 ){
			error("Cannot unpack subsitution model\n");
	}
	else {
		 error("Cannot unpack subsitution model\n");
	}
	
	char *id       = SML_get_data_of_child( node, "id");
	if ( id != NULL ) {
		while( *id < '0' || *id > '9' ){
			id++;
		}
		m->id = atoi( id );
	}
	else m->id = 0;
	
	return m;
}*/


/******************************************************************************************************************************************************/

#pragma mark -
#pragma mark Misc


void update_eigen_system( SubstitutionModel *m ){
	
	make_zero_rows( m->Q, m->nstate);
	normalize_Q( m->Q, m->_freqs, m->nstate );
    
    /*fprintf(stderr, "\n");
    print_frequencies(stderr, m);
    print_rates(stderr, m);
    print_dmatrix(stderr, m->Q, m->nstate, m->nstate, ' ');
    for (int i = 0; i < m->nstate; i++) {
        fprintf(stderr, "%e\n", Parameters_value(m->freqs, i));
    }*/
    
	EigenDecomposition_decompose(m->Q, m->eigendcmp);
    /*printf("Eigen value\n");
    for ( int i = 0; i < m->eigendcmp->dim; i++ ) {
        printf("%d %e\n",i,m->eigendcmp->eval[i]);
    }
    printf("Eigen vector\n");
    for ( int i = 0; i < m->eigendcmp->dim; i++ ) {
        for ( int j = 0; j < m->eigendcmp->dim; j++ ) {
            printf("%e ",m->eigendcmp->evec[i][j]);
        }
        printf("\n\n");
    }*/
    
}

void normalize_Q( double **q, const double *freqs, const int dim ) {
	double subst = 0.0;
	
	for (int i = 0; i < dim; i++) {
		subst += -q[i][i]*freqs[i];
	}
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			q[i][j] /= subst;
		}
	}
}

void make_zero_rows( double **q, const int dim ) {
	for (int i = 0; i < dim; i++) {
		double sum = 0.0;
		for (int j = 0; j < dim; j++) {
			if (i != j) {
				sum += q[i][j];
			}
		}
		q[i][i] = -sum;
	}
}

void check_frequencies( const double *freqs, const int dim ){
	double sum = 0.;
	for ( int i = 0; i < dim; i++)
		sum += freqs[i];
	if( fabs(sum - 1.0) > 0.0001 ){
		fprintf(stderr, "Frequencies do not add up to 1:%f (diff = %f)\n",sum, fabs(sum - 1.0) );
		exit(1);
	}
}

void check_frequencies_parameters( const Parameters* freqs ){
	double sum = 0.;
	for ( int i = 0; i < Parameters_count(freqs); i++)
		sum += Parameters_value(freqs, i);
	if( fabs(sum - 1.0) > 0.0001 ){
		fprintf(stderr, "Frequencies do not add up to 1:%f (diff = %f)\n",sum, fabs(sum - 1.0) );
		exit(1);
	}
}



void bufferize_frequencies(StringBuffer *buffer, const SubstitutionModel *m ){
    if( m->_freqs != NULL ){
		StringBuffer_append_format(buffer, "   A: %f\n", m->_freqs[0] );
		StringBuffer_append_format(buffer, "   C: %f\n", m->_freqs[1] );
		StringBuffer_append_format(buffer, "   G: %f\n", m->_freqs[2] );
		StringBuffer_append_format(buffer, "   T: %f\n", m->_freqs[3] );
    }
    else {
		StringBuffer_append_format(buffer, "   A: 0.25\n" );
		StringBuffer_append_format(buffer, "   C: 0.25\n" );
		StringBuffer_append_format(buffer, "   G: 0.25\n" );
		StringBuffer_append_format(buffer, "   T: 0.25\n" );
    }
}

void bufferize_codon_frequencies(StringBuffer *buffer, const SubstitutionModel *m ){
    char codon[4];
    codon[3] = '\0';
    for ( int i = 0; i < m->nstate; i++ ) {
        GenticCode_encoding_to_codon_string( i, m->gen_code, codon );
        StringBuffer_append_format(buffer, "   %s: %f",  codon, m->_freqs[i]);
        
        if( (i+1) % 8 == 0 ){
            StringBuffer_append_char(buffer, '\n');
        }
    }
}

void bufferize_aa_frequencies(StringBuffer *buffer, const SubstitutionModel *m ){
    char AMINO_ACIDS[20] = "ACDEFGHIKLMNPQRSTVWY";
    
    for ( int i = 0; i < m->nstate; i++ ) {
        StringBuffer_append_format(buffer, "   %c: %f",  AMINO_ACIDS[i], m->_freqs[i]);
        
        if( (i+1) % 5 == 0 ){
            StringBuffer_append_char(buffer, '\n');
        }
    }
}

/*void bufferize_codon_frequencies2(StringBuffer *buffer, const SubstitutionModel *m ){
    char codon[4];
    codon[3] = '\0';
    printf("sdfasD\n");
    for ( int i = 0; i < m->nstate; i++ ) {
        //printf("%d %f\n", i,m->_freqs[i]);
        StringBuffer_append_format(buffer, ",%f", m->_freqs[i]);
        
        
    }
}*/

void print_frequencies(FILE *pfile, const SubstitutionModel *m ){
	StringBuffer *buffer = new_StringBuffer(100);
	if( m->dtype == DATA_TYPE_NUCLEOTIDE ) bufferize_frequencies(buffer, m);
    else if( m->dtype == DATA_TYPE_CODON ) bufferize_codon_frequencies(buffer, m);
    else if( m->dtype == DATA_TYPE_AMINO_ACID ) bufferize_aa_frequencies(buffer, m);
	fprintf(pfile, "%s", buffer->c);
	free_StringBuffer(buffer);
}

/*
 
    A C G T
 A  * a b c
 C  a * d e
 G  b d * 1
 T  c e 1 *
 
 */
void print_rates(FILE *pfile, const SubstitutionModel *m){	
	StringBuffer *buffer = new_StringBuffer(100);
	bufferize_rates(buffer, m);
	fprintf(pfile, "%s", buffer->c);
	free_StringBuffer(buffer);
}

void print_P( double *p, int nstate ){
    for ( int i = 0; i < nstate; i++ ) {
        double tot = 0;
        for ( int j = 0; j < nstate; j++ ) {
            printf("%e ", p[i*nstate+j]);
            tot += p[i*nstate+j];
        }
        printf("| %e\n", tot);
    }
    for ( int i = 0; i < nstate; i++ ) {
        double tot = 0;
        for ( int j = 0; j < nstate; j++ ) {
            tot += p[j*nstate+i];
        }
        printf("%e ", tot);
    }
    printf("\n");
}

void print_P2( double **p, int nstate ){
    for ( int i = 0; i < nstate; i++ ) {
        double tot = 0;
        for ( int j = 0; j < nstate; j++ ) {
            printf("%e ", p[i][j]);
            tot += p[i][j];
        }
        printf("| %e\n", tot);
    }
    for ( int i = 0; i < nstate; i++ ) {
        double tot = 0;
        for ( int j = 0; j < nstate; j++ ) {
            tot += p[j][i];
        }
        printf("%e ", tot);
    }
    printf("\n");
}

void print_substitution_matrix( SubstitutionModel *sm, double t ){
	double **P = dmatrix(sm->nstate, sm->nstate);
	
	fprintf(stderr, "\nBranch length %e\n",t);
	
	int i,j,k,l;
	
	for ( i = 0; i < sm->nstate; i++){
		for ( j = 0; j < sm->nstate; j++){
				for ( k = 0; k < sm->nstate; k++)
					P[i][j] += exp( sm->eigendcmp->eval[k] * t ) * sm->eigendcmp->evec[i][k] * sm->eigendcmp->Invevec[k][j];
		}
	}
//	
//	double temp;
//	for ( i = 0; i < sm->nstate; i++){
//		temp = exp(sm->eigendcmp->eval[i] * t);
//		for ( j = 0; j < sm->nstate; j++){
//			sm->PP[i][j] = sm->eigendcmp->Invevec[i][j] * temp;
//		}
//	}
//	
//	for ( i = 0; i < sm->nstate; i++){
//		for ( j = 0; j < sm->nstate; j++){
//			P[l] = 0.;
//			for ( k = 0; k < sm->eigendcmp->dim; k++ )
//				P[l] += sm->PP[k][j] * sm->eigendcmp->evec[i][k];
//			P[l] = fabs(P[l]);
//			l++;
//		}
//	}
	
	l=0;
	for ( i = 0; i < sm->nstate; i++) {
		for ( j = 0; j < sm->nstate; j++){
			fprintf(stderr, "%e\t", P[i][j]);
			l++;
		}
		fprintf(stderr, "\n");
	}
	
//	fprintf(stderr, "\n");
//	for ( i = 0; i < sm->nstate; i++) {
//		for ( j = 0; j < sm->nstate; j++){
//			fprintf(stderr, "%f\t", sm->Q[i][j]);
//			l++;
//		}
//		fprintf(stderr, "\n");
//	}

	free(P);
}

void bufferize_rates( StringBuffer *buffer, const SubstitutionModel *m ){
	if ( m->dtype == DATA_TYPE_NUCLEOTIDE ) {
		if ( m->modeltype == GTR ) {
			char *rates[] = {
				"AC",
				"AG",
				"AT",
				"CG",
				"CT",
				"GT",
			};
			int i = 0;
			for ( ; i < Parameters_count(m->rates); i++) {
				StringBuffer_append_format(buffer, "   %s: %f\n", rates[i], Parameters_value(m->rates, i));
			}
			StringBuffer_append_string(buffer, "   GT: 1\n");
		}
		else if ( m->modeltype == HKY ) {
			StringBuffer_append_format(buffer, "   Kappa: %f (Ts/Tv: %f)\n", Parameters_value(m->rates, 0), kappa_to_tstv(Parameters_value(m->rates, 0), m->_freqs) );
		}
		else if ( m->modeltype == K80 ) {
			StringBuffer_append_format(buffer, "   Kappa: %f (Ts/Tv: %f)\n", Parameters_value(m->rates, 0), Parameters_value(m->rates, 0) );
		}
        else if( m->modeltype == REVERSIBLE_DNA) {
            char *rates[] = {
				"AC",
				"AG",
				"AT",
				"CG",
				"CT",
				"GT",
			};
            char cars[6] = "abcdef";
			for ( int i = 0; i < 5; i++) {
                char c[2];
                c[0] = cars[m->model[i]];
                c[1] = '\0';
				StringBuffer_append_format(buffer, "   %s: %f", rates[ m->model[i] ], Parameters_value(m->rates, m->model[i]));
				StringBuffer_append_format(buffer, " %s\n", c);
			}
			StringBuffer_append_format(buffer, "   GT: 1\n");
        }
	}
	else if( m->dtype == DATA_TYPE_CODON ){
		if ( m->modeltype == GY94 ) {
			StringBuffer_append_format(buffer, "   Kappa: %f\n", Parameters_value(m->rates, 0));
			StringBuffer_append_format(buffer, "   Omega: %f\n", Parameters_value(m->rates, 1));
		}
	}
}

bool Model_is_valid( const char* model ){
	if (model == NULL ) {
		return false;
	}
	
	for ( int i = 0; i < strlen(model); i++) {
		if ( model[i] < 48 || model[i] > 57) {
			return false;
		}
	}
	return true;
}

int nDNARates( const char *code ){	
	int max = 0;
	for ( int i = 0; i < strlen(code); i++) {
		if ( code[i] >= 48) {
			max = code[i]-48;
		}
	}
	return max+1;
}

/*SubstitutionModel * SubstitutionModel_factory( const char* model_string ){
	if (model_string == NULL ) {
		return NULL;
	}
	
	SubstitutionModel *model = NULL;
	if ( strlen(model_string) == 5 ) {
		if ( !Model_is_valid(model_string) ) {
			error("Model invalid\n");
		}
		else {
			//double freqs[4] = {0.25,0.25,0.25,0.25};
		}
	}
	else if( strcmp(model_string, "GTR") == 0 ){
		model = new_GTR();
	}
	else if( strcmp(model_string, "HKY") == 0 ){
		model = new_HKY();
	}
	else if( strcmp(model_string, "JC69") == 0 ){
		model = new_JC69();
	}
	else if( strcmp(model_string, "K80") == 0 ){
		model = new_K80();
	}
	
	return model;
}

SubstitutionModel * SubstitutionModel_nuc_factory_with_values( modeltype modtype, double *rates, double *freqs ){
	SubstitutionModel *model = NULL;
	switch ( modtype ) {
		case JC69:{
			model = new_JC69();
			break;
        }
		case HKY:{
			model = new_HKY_with_values(freqs, rates[0]);
			break;
        }
		case K80:{
			model = new_K80_with_values(rates[0]);
			break;
        }
		case GTR:{
			model = new_GTR_with_values(freqs, rates);
			break;
        }
		default:
			break;
	}
	return model;
}
modeltype SubstitutionModel_nucleotide_code_to_modeltype( const char *model_string, bool equal_freqs ){
    if( strcmp(model_string, "GY94") == 0 ){
		return GY94;
	}
	else if( strcmp(model_string, "MG94") == 0 ){
		return MG94;
	}
    
    // GTR
    if ( strcasecmp(model_string, "GTR") == 0 || (!equal_freqs && (strcasecmp(model_string, "abcde") == 0 || strcmp(model_string, "01234") == 0 )) ) {
        return GTR;
    }
    // HKY or K80
    else if ( strcasecmp(model_string, "HKY") == 0 || strcasecmp(model_string, "K80") == 0 || strcasecmp(model_string, "abaab") == 0 || strcmp(model_string, "01001") == 0 ) {
        // K80 (if frequencies are equal)
        if ( equal_freqs ) {
            return K80;
        }
        // HKY
        else {
            return HKY;
        }
    }
    // JC69
    else if ( strcasecmp(model_string, "JC69") == 0 || strcasecmp(model_string, "aaaaa") == 0 || strcmp(model_string, "00000") == 0 ) {
        return JC69;
    }
    return REVERSIBLE_DNA;
}

modeltype SubstitutionModel_get_code( const char *model_string ){
    
    
	if( strcasecmp(model_string, "GTR") == 0 ){
		return GTR;
	}
	else if( strcasecmp(model_string, "HKY") == 0 ){
		return HKY;
	}
	else if( strcasecmp(model_string, "JC69") == 0 ){
		return JC69;
	}
	else if( strcasecmp(model_string, "K80") == 0 ){
		return K80;
	}
	else if( strcasecmp(model_string, "GY94") == 0 ){
		return GY94;
	}
	else if( strcasecmp(model_string, "MG94") == 0 ){
		return MG94;
	}
	else if( strcasecmp(model_string, "UREV") == 0 ){
		return NON_REVERSIBLE_DNA;
	}
	else if( strcasecmp(model_string, "NSTAT") == 0 ){
		return NON_STATIONARY_DNA;
	}
	else if( strcasecmp(model_string, "WAG") == 0 ){
		return WAG;
	}
	else if( strcasecmp(model_string, "DAYHOFF") == 0 ){
		return DAYHOFF;
	}
	else if( strcasecmp(model_string, "LG") == 0 ){
		return LG;
	}
	
    if( strlen(model_string) == 5 ){
        return REVERSIBLE_DNA;
	}
    
	fprintf(stderr, "Substitution model not yet implemented (or typpe): %s\n", model_string);
	return INFINITY;
}

char * SubstitutionModel_get_string( modeltype modtype ){
	if( modtype == GTR ){
		return "GTR";
	}
	else if( modtype == HKY ){
		return "HKY";
	}
	else if( modtype == JC69 ){
		return "JC69";
	}
	else if( modtype == K80 ){
		return "K80";
    }
    else if( modtype == GY94 ){
        return "GY94";
    }
    else if( modtype == MG94 ){
        return "MG94";
    }
    else if( modtype == NON_REVERSIBLE_DNA ){
        return "UREV";
    }
    else if( modtype == NON_STATIONARY_DNA ){
        return "NSTAT";
    }
    else if( modtype == WAG ){
        return "WAG";
    }
	else if( modtype == DAYHOFF ){
		return "DAYHOFF";
	}
	else if( modtype == LG ){
		return "LG";
	}
	
	error("Substitution model not recognized\n");

	return NULL;
}*/

/****************************************************************************/


double tstv_to_kappa( double tstv, double *freqs ){
    return (tstv*(freqs[0]+freqs[2])*(freqs[1]+freqs[3]))/((freqs[0]*freqs[2])+(freqs[1]*freqs[3]));
}

double kappa_to_tstv( double kappa, double *freqs ){
    return kappa * ((freqs[0]*freqs[2])+(freqs[1]*freqs[3]))/((freqs[0]+freqs[2])*(freqs[1]+freqs[3]));
}

