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

#include "hashtable.h"

#include "sitepattern.h"

#pragma mark Public function definition


static void _substitution_model_handle_change( Model *self, Model *model, int index ){
	SubstitutionModel* m = (SubstitutionModel*)self->obj;
	m->need_update = true;
	m->dQ_need_update = true;
	self->listeners->fire( self->listeners, self, index );
}

static void _substitution_model_store(Model* self){
	Model** models = (Model**)self->data;
	models[0]->store(models[0]); // simplex freqs
	if(models[1] != NULL) models[1]->store(models[1]); // simplex rates
	if(models[2] != NULL) models[2]->store(models[2]); // discrete parameter
	SubstitutionModel* subst = self->obj;
	if (Parameters_count(subst->rates) > 0) {
		Parameters_store(subst->rates);
	}
}

static void _substitution_model_restore(Model* self){
	Model** models = (Model**)self->data;
	models[0]->restore(models[0]); // simplex freqs
	if(models[1] != NULL) models[1]->restore(models[1]); // simplex rates
	if(models[2] != NULL) models[2]->restore(models[2]); // discrete parameter
	SubstitutionModel* subst = self->obj;
	if (Parameters_count(subst->rates) > 0) {
		bool changed = false;
		Parameter*p = NULL;
		for (int i = 0; i < Parameters_count(subst->rates); i++) {
			p = Parameters_at(subst->rates, i);
			if (Parameter_changed(p)) {
				changed = true;
			}
			Parameter_restore_quietly(p);
		}
		if (changed) {
			p->listeners->fire_restore(p->listeners, NULL, p->id);
		}
	}
}

static void _substitution_model_handle_restore( Model *self, Model *model, int index ){
	SubstitutionModel* m = (SubstitutionModel*)self->obj;
	m->need_update = true;
	m->dQ_need_update = true;
	self->listeners->fire_restore( self->listeners, self, index );
}

static void _substitution_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free subsitution model %s\n", self->name);
		SubstitutionModel* m = (SubstitutionModel*)self->obj;
		if(self->data != NULL){
			Model** models = (Model**)self->data;
			for (int i = 0; i < 3; i++) {
				if(models[i] != NULL){
					models[i]->free(models[i]);
				}
			}
			free(models);
		}
		
		free_DataType(m->datatype);
		free(m->name);
		if( m->eigendcmp != NULL ) free_EigenDecomposition(m->eigendcmp);
		if( m->Q != NULL ) free_dmatrix(m->Q, m->nstate);
		if( m->PP != NULL ) free_dmatrix(m->PP, m->nstate);
		if( m->rates != NULL ) free_Parameters(m->rates);
		if( m->model != NULL ) m->model->free(m->model);
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
	
	Model* msimplex_freq_clone = NULL;
	Model* msimplex_rates_clone = NULL;
	
	if(self->data != NULL){
		Model** msimplex = (Model**)self->data;
		if(subst->simplex != NULL){
			if (Hashtable_exists(hash, msimplex[0]->name)) {
				msimplex_freq_clone = Hashtable_get(hash, msimplex[0]->name);
				msimplex_freq_clone->ref_count++;
			}
			else{
				msimplex_freq_clone = msimplex[0]->clone(msimplex[0], hash);
				Hashtable_add(hash, msimplex_freq_clone->name, msimplex_freq_clone);
			}
		}
		
		if(subst->rates_simplex != NULL){
			if (Hashtable_exists(hash, msimplex[1]->name)) {
				msimplex_rates_clone = Hashtable_get(hash, msimplex[1]->name);
				msimplex_rates_clone->ref_count++;
			}
			else{
				msimplex_rates_clone = msimplex[1]->clone(msimplex[1], hash);
				Hashtable_add(hash, msimplex_rates_clone->name, msimplex_rates_clone);
			}
		}
		//TODO: discrete parameter
	}
	
	Parameters* ps = NULL;
	Simplex* rates_simplex = NULL;
	if(msimplex_rates_clone == NULL && Parameters_count(subst->rates) > 0){
		ps = new_Parameters(1);
		for (int i = 0; i < Parameters_count(subst->rates); i++) {
			char* name = Parameters_name(subst->rates, i);
			if (Hashtable_exists(hash, name)) {
				Parameters_add(ps, Hashtable_get(hash, name));
			}
			else{
				Parameter* p = clone_Parameter(Parameters_at(subst->rates, i));
				Parameters_move(ps, p);
				Hashtable_add(hash, name, p);
			}
		}
	}
	SubstitutionModel* mclone = NULL;
	if(msimplex_freq_clone != NULL){
		mclone = clone_substitution_model_with(subst, ps, (Simplex*)msimplex_freq_clone->obj);
	}
	else{
		mclone = clone_substitution_model_with(subst, ps, NULL);
	}
	
	if (msimplex_rates_clone != NULL) {
		mclone->rates_simplex = (Simplex*)msimplex_rates_clone->obj;
	}

	if (mclone->model != NULL) {
		//TODO: use discreteparameter
//		Hashtable_add(hash, mclone->model->name, mclone->model);
	}
	mclone->datatype = subst->datatype;
	mclone->datatype->ref_count++;
	Model *clone = new_SubstitutionModel2(self->name, mclone, msimplex_freq_clone, msimplex_rates_clone);
	Hashtable_add(hash, clone->name, clone);
	if(msimplex_freq_clone != NULL) msimplex_freq_clone->free(msimplex_freq_clone);
	free_Parameters(ps);
	return clone;
}

Model * new_SubstitutionModel3( const char* name, SubstitutionModel *sm, Model* freqs_simplex, Model* rates_simplex, Model* discrete_model ){
	Model *model = new_Model(MODEL_SUBSTITUTION, name, sm);
	Model** models = (Model**)calloc(3, sizeof(Model*));
	for(size_t i = 0; i < 3; i++){
		models[i] = NULL;
	}
	model->data = models;
	if(freqs_simplex != NULL){
		models[0] = freqs_simplex;
		freqs_simplex->ref_count++;
		freqs_simplex->listeners->add( freqs_simplex->listeners, model );
	}
	
	int i = 0;
	if ( sm->rates != NULL ) {
		for ( i = 0; i < Parameters_count(sm->rates); i++ ) {
			Parameters_at(sm->rates, i)->listeners->add( Parameters_at(sm->rates, i)->listeners, model );
		}
	}
	
	if (rates_simplex != NULL) {
		models[1] = rates_simplex;
		rates_simplex->ref_count++;
		rates_simplex->listeners->add( rates_simplex->listeners, model );
	}
	
	if (discrete_model != NULL) {
		models[2] = discrete_model;
		discrete_model->ref_count++;
		discrete_model->listeners->add( discrete_model->listeners, model );
	}
	
	model->update = _substitution_model_handle_change;
	model->handle_restore = _substitution_model_handle_restore;
	model->store = _substitution_model_store;
	model->restore = _substitution_model_restore;
	model->free = _substitution_model_free;
	model->clone = _substitution_model_clone;
	return model;
}


// SubstitutionModel2 listen to the rate and freq parameters
Model * new_SubstitutionModel2( const char* name, SubstitutionModel *sm, Model* freqs_simplex, Model* rates_simplex ){
	return new_SubstitutionModel3(name, sm, freqs_simplex, rates_simplex, NULL);
}

Model* new_SubstitutionModel_from_json(json_node* node, Hashtable*hash){
	char* allowed[] ={
		"datatype",
		"frequencies",
		"init",
		"model",
		"normalize",
		"rates",
		"structure"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	json_node* model_node = get_json_node(node, "model");
	json_node* freqs_node = get_json_node(node, "frequencies");
	json_node* rates_node = get_json_node(node, "rates");
	json_node* datatype_node = get_json_node(node, "datatype");
	json_node* structure_node = get_json_node(node, "structure");
	json_node* init_node = get_json_node(node, "init");
	bool normalize = get_json_node_value_bool(node, "normalize", true);
	char* id = get_json_node_value_string(node, "id");
	
	// DataType
	DataType* datatype = NULL;
	
	if (datatype_node->node_type == MJSON_STRING && ((char*)datatype_node->value)[0] == '&') {
		char* ref = (char*)datatype_node->value;
		// check it starts with a &
		datatype = Hashtable_get(hash, ref+1);
		datatype->ref_count++;
	}
	else{
		datatype = new_DataType_from_json(datatype_node, hash);
	}

	// Simplex Frequencies
	Model* mfreqs_simplex = NULL;
	Simplex* freqs_simplex = NULL;
	if (freqs_node!= NULL && freqs_node->node_type == MJSON_OBJECT) {
		mfreqs_simplex = new_SimplexModel_from_json(freqs_node, hash);
		json_node* id_node = get_json_node(freqs_node, "id");
		Hashtable_add(hash, (char*)id_node->value, mfreqs_simplex);
		freqs_simplex =  (Simplex*)mfreqs_simplex->obj;
		for (int i = 0; i < Parameters_count(freqs_simplex->parameters); i++) {
			Parameters_at(freqs_simplex->parameters, i)->model = MODEL_SUBSTITUTION;
		}
	}
	else if(freqs_node!= NULL && freqs_node->node_type == MJSON_STRING){
		char* ref = (char*)freqs_node->value;
		// check it starts with a &
		mfreqs_simplex = Hashtable_get(hash, ref+1);
		mfreqs_simplex->ref_count++;
		freqs_simplex =  (Simplex*)mfreqs_simplex->obj;
	}
	
	// Rates
	Parameters* rates = NULL;
	Model* mrates_simplex = NULL;
	Simplex* rates_simplex = NULL;
	
	// Rates is made of doubles
	if (rates_node != NULL && rates_node->node_type == MJSON_ARRAY) {
		rates = new_Parameters_from_json(rates_node, hash);
	}
	else if (rates_node != NULL && rates_node->node_type == MJSON_OBJECT) {
		char* rate_simplex_ref = get_json_node_value_string(rates_node, "type");
		// Rates is a simplex
		if (rate_simplex_ref != NULL) {
			if (rates_node->node_type == MJSON_OBJECT) {
				mrates_simplex = new_SimplexModel_from_json(rates_node, hash);
				char* id = get_json_node_value_string(rates_node, "id");
				Hashtable_add(hash, id, mrates_simplex);
			}
			else if(rates_node->node_type == MJSON_STRING){
				char* ref = (char*)rates_node->value;
				// check it starts with a &
				mrates_simplex = Hashtable_get(hash, ref+1);
				mrates_simplex->ref_count++;
			}
			rates_simplex = (Simplex*)mrates_simplex->obj;
			for (int i = 0; i < Parameters_count(rates_simplex->parameters); i++) {
				Parameters_at(rates_simplex->parameters, i)->model = MODEL_SUBSTITUTION;
			}
		}
		// Rates is made of parameters
		else{
			rates = new_Parameters(rates_node->child_count);
			for(int i = 0; i < rates_node->child_count; i++){
				json_node* p_node = rates_node->children[i];
				Parameter* p = new_Parameter_from_json(p_node, hash);
				p->model = MODEL_SUBSTITUTION;
				Parameters_move(rates, p);
			}
		}
	}
	for (int i = 0; i < Parameters_count(rates); i++) {
		Hashtable_add(hash, Parameters_name(rates, i), Parameters_at(rates, i));
	}
	
	// Rate assignment
	char** assignment = NULL;
	if (Parameters_count(rates) > 0) {
		assignment = (char**)malloc(sizeof(char*)*rates_node->child_count);
		for (int i = 0; i < rates_node->child_count; i++) {
			json_node* child = rates_node->children[i];
			assignment[i] = child->key;
		}
	}
	
	// General model
	//TODO: use rate simplex
	if (structure_node != NULL) {
		Model* mdp = new_DiscreteParameterModel_from_json(model_node, hash);
		check_constraints(rates, 0, INFINITY, 0.001, 1000);
		SubstitutionModel* m = new_GeneralModel_with_parameters(datatype, mdp->obj, rates, freqs_simplex, -1, normalize);
		Model* mm = new_SubstitutionModel3(id, m, mfreqs_simplex, mrates_simplex, mdp);

		free_DataType(datatype);
		free_Parameters(rates);
		mdp->free(mdp);
		if(mfreqs_simplex != NULL) mfreqs_simplex->free(mfreqs_simplex);
		if(mrates_simplex != NULL) mrates_simplex->free(mrates_simplex);
		if(assignment != NULL) free(assignment);
		return mm;
	}
	
	char* model_string = (char*)model_node->value;
	
	if (datatype->type == DATA_TYPE_GENERIC) {
		// datatype->ref_count already incremented
		SubstitutionModel *m = new_GeneralJC69Model_with_parameters( datatype, freqs_simplex, normalize );
		Model* mm = new_SubstitutionModel2(id, m, mfreqs_simplex, mrates_simplex);
		if(mfreqs_simplex != NULL) mfreqs_simplex->free(mfreqs_simplex);
		return mm;
	}
	
	SubstitutionModel* m = SubstitutionModel_factory(model_string, datatype, freqs_simplex, rates_simplex, rates, (const char**)assignment);
	
	if(mfreqs_simplex == NULL  && m->simplex != NULL){
		mfreqs_simplex = new_SimplexModel("anonymousfreqs", m->simplex);
		for (int i = 0; i < Parameters_count(m->simplex->parameters); i++) {
			Parameters_at(m->simplex->parameters, i)->model = MODEL_SUBSTITUTION;
		}
	}
	if(init_node != NULL && datatype->type == DATA_TYPE_NUCLEOTIDE && (Parameters_count(m->rates) >= 5 || m->rates_simplex != NULL)){
		json_node* patterns_node = get_json_node(init_node, "sitepattern");
		char* patterns_ref = (char*)patterns_node->value;
		SitePattern* patterns = Hashtable_get(hash, patterns_ref+1);
		double** mat = SitePattern_rates(patterns);
		double* v = dvector(6);
		int index = 0;
		for (int i = 0; i < 4; i++) {
			for (int j = i+1; j < 4; j++) {
				v[index++] = mat[i][j] + mat[j][i];
			}
		}
		if ( m->relativeTo == -1) {
			double sum = 0;
			for (int i = 0; i < 6; i++) {
				sum += v[i];
			}
			for (int i = 0; i < 6; i++) {
				v[i] /= sum;
			}
		}
		else{
			for (int i = 0; i < 6; i++) {
				v[i] /= v[m->relativeTo];
			}
			int moveCount = 6 - m->relativeTo - 1;
			if(moveCount > 0){
				memcpy(v+m->relativeTo, v+m->relativeTo+1, moveCount*sizeof(double));
			}
		}
		if(m->rates_simplex != NULL){
			m->rates_simplex->set_values(m->rates_simplex, v);
		}
		else{
			m->set_rates(m, v);
		}
		for (int i = 0; i < 4; i++) {
			v[i] = mat[i][i];
		}
		m->simplex->set_values(m->simplex, v);
		free(v);
		free_dmatrix(mat, datatype->state_count(datatype));
	}
	
	free_DataType(datatype);
	Model* mm = new_SubstitutionModel2(id, m, mfreqs_simplex, mrates_simplex);
	free_Parameters(rates);
	if(mfreqs_simplex != NULL) mfreqs_simplex->free(mfreqs_simplex);
	if(mrates_simplex != NULL) mrates_simplex->free(mrates_simplex);
	if(assignment != NULL) free(assignment);
	return mm;
}

double get_frequency( SubstitutionModel *m, int base ){
	return m->simplex->get_value(m->simplex, base);
}

const double * get_frequencies( SubstitutionModel *m ){
	return m->simplex->get_values(m->simplex);
}


// not normalized
void build_Q_flat(double* Q, const double* freqs, size_t stateCount){
    for ( size_t i = 0; i < stateCount; i++ )  {
        for ( size_t j = 0; j < stateCount; j++ ) {
            Q[i*stateCount+j] *= freqs[j];
        }
    }
    
    for ( size_t i = 0; i < stateCount; i++ )  {
		double sum = 0.0;
		Q[i*stateCount+i] = 0;
		for ( size_t j = 0; j < stateCount; j++ ) {
			sum += Q[i*stateCount+j];
		}
		Q[i*stateCount+i] = -sum;
	}
}

//MARK: derivatives

void dPdp_with_dQdp(SubstitutionModel *m, double* dQ, double* dP, double t){
	size_t stateCount = m->nstate;
	double *v = m->eigendcmp->eval;
	double* xx = dvector(stateCount*stateCount);
	Matrix_mult3(dP, (const double**)m->eigendcmp->Invevec, m->dQ, stateCount, stateCount, stateCount, stateCount);
	Matrix_mult4(xx, dP, (const double**)m->eigendcmp->evec, stateCount, stateCount, stateCount, stateCount);
	// up to now the above operations can be recycled across branches
	for(size_t i = 0; i < stateCount; i++){
		for(size_t j = 0; j < stateCount; j++){
			if(v[i] != v[j]){
				dP[i*stateCount+j] = xx[i*stateCount+j]*(exp(v[i]*t) - exp(v[j]*t))/(v[i]-v[j]);
			}
			else{
				dP[i*stateCount+j] = xx[i*stateCount+j]*t*exp(v[i]*t);
			}
		}
	}
	Matrix_mult3(xx, (const double**)m->eigendcmp->evec, dP, stateCount, stateCount, stateCount, stateCount);
	Matrix_mult4(dP, xx, (const double**)m->eigendcmp->Invevec, stateCount, stateCount, stateCount, stateCount);
	free(xx);
}

#pragma mark -
#pragma mark Private functions

static void foo_update( SubstitutionModel *m ){} // does not nothing

static void _set_rates( SubstitutionModel *m, const double *rates ){
	for ( int i = 0; i < Parameters_count(m->rates); i++) {
		Parameters_set_value(m->rates,i, rates[i]);
	}
	m->need_update = true;
    m->dQ_need_update = true;
}

// Set the relative frequency (ie. Parameter)!!
static void _set_relative_frequency( SubstitutionModel *m, const double freq, const int index ){
	m->simplex->set_parameter_value(m->simplex, index, freq);
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
		update_eigen_system(m);
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
		update_eigen_system(m);
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
		update_eigen_system(m);
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
		update_eigen_system(m);
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
		update_eigen_system(m);
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
		update_eigen_system(m);
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
		update_eigen_system(m);
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
		update_eigen_system(m);
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
		update_eigen_system(m);
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

#pragma mark -


SubstitutionModel * create_substitution_model( const char *name, const modeltype modelname, DataType* datatype, Simplex* freqs ){
    SubstitutionModel *m = (SubstitutionModel *)malloc( sizeof(SubstitutionModel) );
    assert(m);
    
    m->nstate = 0;
    m->name = String_clone(name);
    m->modeltype = modelname;
    m->datatype = datatype;
    
    // Parameters and Q matrix
    m->Q = NULL;
    m->PP = NULL;
    
    m->rates = NULL;
	m->rates_simplex = NULL;
	m->simplex = freqs;
	
    m->eigendcmp = NULL;
    
    m->need_update = true;
    m->dQ_need_update = true;
    
	// Functions
	m->get_frequencies = get_frequencies;
	
	m->set_rates = _set_rates;
	
	m->set_relative_frequency = _set_relative_frequency; // relative frequency with dimention nstate-1
    m->set_rate = _set_rate;
    m->update_Q = foo_update;
    
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
	m->relativeTo = -1;
    m->grad_wrt_reparam = true;
    
    return m;
}

SubstitutionModel * create_nucleotide_model( const char *name, const modeltype modelname, Simplex* freqs ){
	DataType* datatype = new_NucleotideDataType();
    SubstitutionModel *m = create_substitution_model(name, modelname, datatype, freqs);
    m->nstate = 4;
    m->Q = dmatrix(m->nstate, m->nstate);
    m->PP = dmatrix(m->nstate, m->nstate);
    m->eigendcmp = new_EigenDecomposition(m->nstate);
    return m;
}

SubstitutionModel * create_codon_model( const char *name, const modeltype modelname, unsigned gen_code, Simplex* freqs ){
	DataType* datatype = new_CodonDataType(gen_code);
    SubstitutionModel *m = create_substitution_model(name, modelname, datatype, freqs);
    m->nstate = NUMBER_OF_CODONS[gen_code];
    m->gen_code = gen_code;
    m->Q = dmatrix(m->nstate, m->nstate);
    m->PP = dmatrix(m->nstate, m->nstate);
    m->eigendcmp = new_EigenDecomposition(m->nstate);
    return m;
}

SubstitutionModel * create_aa_model( const char *name, const modeltype modelname, Simplex* freqs ){
	DataType* datatype = new_AminoAcidDataType();
    SubstitutionModel *m = create_substitution_model(name, modelname, datatype, freqs);
    m->nstate = 20;
    m->Q = dmatrix(m->nstate, m->nstate);
    m->PP = dmatrix(m->nstate, m->nstate);
    m->eigendcmp = new_EigenDecomposition(m->nstate);
    return m;
}

SubstitutionModel * create_general_model( const char *name, const modeltype modelname, DataType* datatype, Simplex* freqs){
    SubstitutionModel *m = create_substitution_model(name, modelname, datatype, freqs);
    m->nstate = freqs->K;
    m->Q = dmatrix(m->nstate, m->nstate);
    m->PP = dmatrix(m->nstate, m->nstate);
    m->eigendcmp = new_EigenDecomposition(m->nstate);
    return m;
}

void free_SubstitutionModel( SubstitutionModel *m){
	free(m->name);
	if( m->eigendcmp != NULL ) free_EigenDecomposition(m->eigendcmp);
	if( m->Q != NULL ) free_dmatrix(m->Q, m->nstate);
	if( m->PP != NULL ) free_dmatrix(m->PP, m->nstate);
	if( m->rates != NULL ) free_Parameters(m->rates);
    if( m->model != NULL ) m->model->free(m->model);
    if( m->dQ != NULL ) free(m->dQ);
	if(m->simplex != NULL) free_Simplex(m->simplex);
	if(m->rates_simplex != NULL) free_Simplex(m->rates_simplex);
	free_DataType(m->datatype);
	free(m);
}


SubstitutionModel * clone_substitution_model(SubstitutionModel *m){
	SubstitutionModel *clone =  create_substitution_model( m->name, m->modeltype, m->datatype, clone_Simplex(m->simplex));
	clone->nstate = m->nstate;
	if (m->rates_simplex != NULL) {
		clone->rates_simplex = clone_Simplex(m->rates_simplex);
	}
	if( m->Q != NULL ){
		clone->eigendcmp = clone_EigenDecomposition(m->eigendcmp);
		clone->Q = clone_dmatrix( m->Q, m->nstate, m->nstate );
		clone->PP = clone_dmatrix( m->PP, m->nstate, m->nstate );
	}
	
	clone->gen_code = m->gen_code;
	clone->reversible = m->reversible;
	clone->normalize = m->normalize;
	
	// Functions
	clone->p_t   = m->p_t;
	clone->p_t_transpose   = m->p_t_transpose;
	
	clone->dp_dt   = m->dp_dt;
	clone->dp_dt_transpose   = m->dp_dt_transpose;
	
	clone->d2p_d2t   = m->d2p_d2t;
	clone->d2p_d2t_transpose   = m->d2p_d2t_transpose;
	
	clone->dPdp = m->dPdp;
	if( m->dQ != NULL ) clone->dQ = clone_dvector(m->dQ, m->nstate*m->nstate);
	
	clone->update_Q = m->update_Q;
	clone->set_rates = m->set_rates;
	clone->set_relative_frequency = m->set_relative_frequency;
	clone->set_rate = m->set_rate;
	
	clone->free = m->free;
	clone->clone = m->clone;
	
	if( m->rates != NULL ){
		clone->rates = clone_Parameters(m->rates );
		
	}
	
	if( m->model != NULL ) clone->model = m->model->clone(m->model);
	
	clone->need_update = false;
	clone->dQ_need_update = false;
	clone->relativeTo = m->relativeTo;

	clone->grad_wrt_reparam = m->grad_wrt_reparam;
	
	return clone;
}

SubstitutionModel * clone_substitution_model_with(SubstitutionModel *m, const Parameters* rates, Simplex* simplex){
	SubstitutionModel *clone =  create_substitution_model( m->name, m->modeltype, m->datatype, simplex);
	clone->nstate = m->nstate;
	clone->simplex = simplex;
	
	if(rates != NULL){
		clone->rates = new_Parameters(Parameters_count(rates));
		Parameters_add_parameters(clone->rates, rates);
	}
	
	if( m->Q != NULL ){
		clone->eigendcmp = clone_EigenDecomposition(m->eigendcmp);
		clone->Q = clone_dmatrix( m->Q, m->nstate, m->nstate );
		clone->PP = clone_dmatrix( m->PP, m->nstate, m->nstate );
	}
	
	clone->gen_code = m->gen_code;
	clone->reversible = m->reversible;
	clone->normalize = m->normalize;
	
	// Functions
	clone->p_t   = m->p_t;
	clone->p_t_transpose   = m->p_t_transpose;
	
	clone->dp_dt   = m->dp_dt;
	clone->dp_dt_transpose   = m->dp_dt_transpose;
	
	clone->d2p_d2t   = m->d2p_d2t;
	clone->d2p_d2t_transpose   = m->d2p_d2t_transpose;
	
	clone->dPdp = m->dPdp;
	if( m->dQ != NULL ) clone->dQ = clone_dvector(m->dQ, m->nstate*m->nstate);
	
	clone->update_Q = m->update_Q;
	clone->set_rates = m->set_rates;
	clone->set_relative_frequency = m->set_relative_frequency;
	clone->set_rate = m->set_rate;
	
	clone->free = m->free;
	clone->clone = m->clone;
	
	if( m->model != NULL ) clone->model = m->model->clone(m->model);
	
	clone->need_update = true;
	clone->dQ_need_update = true;
	clone->relativeTo = m->relativeTo;

	clone->grad_wrt_reparam = m->grad_wrt_reparam;
	
	return clone;
}

/******************************************************************************************************************************************************/

#pragma mark -
#pragma mark Misc


void update_eigen_system( SubstitutionModel *m ){
    
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

double normalizing_constant_Q( double **q, const double *freqs, size_t dim ) {
	double subst = 0.0;
	
	for (size_t i = 0; i < dim; i++) {
		subst += -q[i][i]*freqs[i];
	}
	return subst;
}

double normalizing_constant_Q_flat( double *q, const double *freqs, size_t dim ) {
	double subst = 0.0;
	
	for (size_t i = 0; i < dim; i++) {
		subst += -q[i*dim+i]*freqs[i];
	}
	return subst;
}

double normalize_Q( double **q, const double *freqs, size_t dim ) {
	double normalizingConstant = normalizing_constant_Q(q, freqs, dim);
	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < dim; j++) {
			q[i][j] /= normalizingConstant;
		}
	}
	return normalizingConstant;
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

void check_frequencies_parameters( Parameters* freqs ){
	double sum = 0.;
	for ( int i = 0; i < Parameters_count(freqs); i++)
		sum += Parameters_value(freqs, i);
	if( fabs(sum - 1.0) > 0.0001 ){
		fprintf(stderr, "Frequencies do not add up to 1:%f (diff = %f)\n",sum, fabs(sum - 1.0) );
		exit(1);
	}
}



void bufferize_frequencies(StringBuffer *buffer, SubstitutionModel *m ){
	const double* freqs = m->get_frequencies(m);
	StringBuffer_append_format(buffer, "   A: %f\n", freqs[0] );
	StringBuffer_append_format(buffer, "   C: %f\n", freqs[1] );
	StringBuffer_append_format(buffer, "   G: %f\n", freqs[2] );
	StringBuffer_append_format(buffer, "   T: %f\n", freqs[3] );
}

void bufferize_codon_frequencies(StringBuffer *buffer, SubstitutionModel *m ){
    char codon[4];
    codon[3] = '\0';
	const double* freqs = m->get_frequencies(m);
    for ( int i = 0; i < m->nstate; i++ ) {
        GenticCode_encoding_to_codon_string( i, m->gen_code, codon );
        StringBuffer_append_format(buffer, "   %s: %f",  codon, freqs[i]);
        
        if( (i+1) % 8 == 0 ){
            StringBuffer_append_char(buffer, '\n');
        }
    }
}

void bufferize_aa_frequencies(StringBuffer *buffer, SubstitutionModel *m ){
    char AMINO_ACIDS[20] = "ACDEFGHIKLMNPQRSTVWY";
    const double* freqs = m->get_frequencies(m);
    for ( int i = 0; i < m->nstate; i++ ) {
        StringBuffer_append_format(buffer, "   %c: %f",  AMINO_ACIDS[i], freqs[i]);
        
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

void print_frequencies(FILE *pfile, SubstitutionModel *m ){
	StringBuffer *buffer = new_StringBuffer(100);
	if( m->datatype->type == DATA_TYPE_NUCLEOTIDE ) bufferize_frequencies(buffer, m);
    else if( m->datatype->type == DATA_TYPE_CODON ) bufferize_codon_frequencies(buffer, m);
    else if( m->datatype->type == DATA_TYPE_AMINO_ACID ) bufferize_aa_frequencies(buffer, m);
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
	if ( m->datatype->type == DATA_TYPE_NUCLEOTIDE ) {
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
			StringBuffer_append_format(buffer, "   Kappa: %f (Ts/Tv: %f)\n", Parameters_value(m->rates, 0), kappa_to_tstv(Parameters_value(m->rates, 0), m->get_frequencies(m)) );
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
                c[0] = cars[m->model->values[i]];
                c[1] = '\0';
				StringBuffer_append_format(buffer, "   %s: %f", rates[i], Parameters_value(m->rates, m->model->values[i]));
				StringBuffer_append_format(buffer, " %s\n", c);
			}
			StringBuffer_append_format(buffer, "   GT: 1\n");
        }
	}
	else if( m->datatype->type == DATA_TYPE_CODON ){
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

SubstitutionModel * SubstitutionModel_factory( const char* model_string, DataType* datatype, Simplex* freqSimplex, Simplex* rates_simplex, const Parameters* rates, const char** assignment ){
	SubstitutionModel *mod = NULL;
	size_t K = datatype->state_count(datatype);
	
	if (datatype->type == DATA_TYPE_NUCLEOTIDE) {
		check_constraints((Parameters*)rates, 0, INFINITY, 0.001, 100);
		bool equal_frequencies = false;
		for (int i = 0; i < K; i++) {
			Parameter* p = Parameters_at(freqSimplex->parameters, i);
			if (Parameter_estimate(p) != false || Parameter_value(p) != 0.25) {
				equal_frequencies = true;
				break;
			}
		}
		char* model_string2 = NULL;
		if( model_string[0] == '0' ){
			if(strcasecmp("00000", model_string) == 0 && equal_frequencies){
				model_string2 = String_clone("JC69");
			}
			else if( strcasecmp("01001", model_string) == 0 && equal_frequencies){
				model_string2 = String_clone("K80");
			}
			else if(strcasecmp("00000", model_string) == 0){
				model_string2 = String_clone("F81");
			}
			else if( strcasecmp("01001", model_string) == 0){
				model_string2 = String_clone("HKY");
			}
			else if( strcasecmp("01234", model_string) == 0){
				model_string2 = String_clone("GTR");
			}
			else{
				model_string2 = String_clone(model_string);
			}
		}
		else{
			model_string2 = String_clone(model_string);
		}
		
		if( strcasecmp("JC69", model_string2) == 0 ){
			mod = new_JC69(freqSimplex);
		}
		else if( strcasecmp("K80", model_string2) == 0 ){
			mod = new_K80_with_parameters(Parameters_at(rates, 0));
		}
		else if( strcasecmp("F81", model_string2) == 0 ){
			mod = new_F81(freqSimplex);
		}
		else if( strcasecmp("HKY", model_string2) == 0 ){
			mod = new_HKY_with_parameters(freqSimplex, Parameters_at(rates, 0));
		}
		else if( strcasecmp("GTR", model_string2) == 0 || strcasecmp("SYM", model_string2) == 0 ){
			if(rates_simplex == NULL){
				Parameter* ac = NULL;
				Parameter* ag = NULL;
				Parameter* at = NULL;
				Parameter* cg = NULL;
				Parameter* ct = NULL;
				Parameter* gt = NULL;
				for (int i = 0; i < Parameters_count(rates); i++) {
					if (strcasecmp(assignment[i], "ac") == 0) {
						ac = Parameters_at(rates, i);
					}
					else if (strcasecmp(assignment[i], "ag") == 0) {
						ag = Parameters_at(rates, i);
					}
					else if (strcasecmp(assignment[i], "at") == 0) {
						at = Parameters_at(rates, i);
					}
					else if (strcasecmp(assignment[i], "cg") == 0) {
						cg = Parameters_at(rates, i);
					}
					else if (strcasecmp(assignment[i], "ct") == 0) {
						ct = Parameters_at(rates, i);
					}
					else if (strcasecmp(assignment[i], "gt") == 0) {
						gt = Parameters_at(rates, i);
					}
				}
				mod = new_GTR_with_parameters(freqSimplex, ac, ag, at, cg, ct, gt);
			}
			else{
				mod= new_GTR_with_simplexes(freqSimplex, rates_simplex);
			}
		}
		else if( strcasecmp("UREV", model_string2) == 0 ){
			mod = new_UnrestrictedNucleotideModel_with_parameters(freqSimplex, rates);
		}
		else if( strcasecmp("NONSTAT", model_string2) == 0 ){
			mod = new_NONSTATNucleotideModel_with_parameters(freqSimplex, rates);
		}
		else{
			mod = new_ReversibleNucleotideModel_with_parameters(model_string2, freqSimplex, rates);
		}
		free(model_string2);
	}
	else if (datatype->type == DATA_TYPE_AMINO_ACID) {
		if ( strcasecmp("LG", model_string) == 0) {
			mod = new_LG_with_parameters(freqSimplex);
		}
		else if ( strcasecmp("WAG", model_string) == 0) {
			mod = new_WAG_with_parameters(freqSimplex);
		}
		else if ( strcasecmp("DAYHOFF", model_string) == 0) {
			mod = new_DAYHOFF_with_parameters(freqSimplex);
		}
		else{
			fprintf(stderr, "Amino acid Subsitution model unknown: %s\n", model_string);
			exit(1);
		}
	}
	else if (datatype->type == DATA_TYPE_CODON) {
		if ( strcasecmp("GY94", model_string) == 0) {
			
		}
		else if ( strcasecmp("MG94", model_string) == 0) {
			
		}
		else{
			fprintf(stderr, "Codon Subsitution model unknown: %s\n", model_string);
			exit(1);
		}
	}
	else if (datatype->type == DATA_TYPE_GENERIC) {
		exit(1);
	}
	else{
		error("factory Susbtition model\n");
	}
	mod->grad_wrt_reparam = true;
	return mod;
}

/****************************************************************************/


double tstv_to_kappa( double tstv, const double *freqs ){
    return (tstv*(freqs[0]+freqs[2])*(freqs[1]+freqs[3]))/((freqs[0]*freqs[2])+(freqs[1]*freqs[3]));
}

double kappa_to_tstv( double kappa, const double *freqs ){
    return kappa * ((freqs[0]*freqs[2])+(freqs[1]*freqs[3]))/((freqs[0]+freqs[2])*(freqs[1]+freqs[3]));
}

