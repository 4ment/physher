/*
 *  substmodel.h
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

#ifndef _SUBSTITUTION_MODEL_H_
#define _SUBSTITUTION_MODEL_H_

#include "sequence.h"
#include "eigen.h"
#include "parameters.h"
#include "simplex.h"

#include "substmodels.h"

typedef enum modeltype {
    // Nucleotide models
	GTR  = 1,
	HKY  = 2,
	JC69 = 3,
	K80  = 4,
	NON_REVERSIBLE_DNA = 5,
	NON_STATIONARY_DNA = 6,
	REVERSIBLE_DNA = 0,
	
    // Codon models
	MG94        = 101,
	GY94        = 102,
	OTHER_CODON = 100,
	
    // Amino acid models
	WAG         = 201,
	FLU         = 202,
	DAYHOFF     = 203,
	LG          = 204,
	OTHER_AA    = 200,
    
    // Generic models
    REVERSIBLE     = 300,
    NONREVERSIBLE  = 400
} modeltype;


typedef struct SubstitutionModel{
	//double *_freqs; // hold the real frequencies ex: 0.25 0.25 0.25 0.25
	Simplex* simplex;
	unsigned int nstate;
	char *name;
	modeltype modeltype;
	datatype dtype;
	EigenDecomposition *eigendcmp;
	double **Q;
	DiscreteParameter* model;
    
    double *dQ;
    bool dQ_need_update;
	
	bool need_update;
	
	Parameters *rates;
	//Parameters *freqs; // hold relative frequencies

	int relativeTo;
	
	// temp variables
	double **PP;
	
	unsigned gen_code;
    bool reversible;
    bool normalize;
	
	void (*update_Q)( struct SubstitutionModel * );
	
//	void (*update_frequencies)( struct SubstitutionModel * );
	
	const double* (*get_frequencies)( struct SubstitutionModel * );
	
	void (*set_rates)( struct SubstitutionModel *, const double * );
	
	void (*set_relative_frequency)( struct SubstitutionModel *, const double, const int );
	void (*set_rate)( struct SubstitutionModel *, const double, const int );
	
	double (*pij_t)( struct SubstitutionModel *, const int, const int, const double );
	void (*p_t)( struct SubstitutionModel *, const double, double * );
	void (*p_t_transpose)( struct SubstitutionModel *, const double, double * );
    
    // Derivatives
	void (*dp_dt)( struct SubstitutionModel *, const double, double * );
	void (*dp_dt_transpose)( struct SubstitutionModel *, const double, double * );
    
	void (*d2p_d2t)( struct SubstitutionModel *, const double, double * );
	void (*d2p_d2t_transpose)( struct SubstitutionModel *, const double, double * );
    
    void (*dPdp)(struct SubstitutionModel *m, int index, double* mat, double t);//derivative of P with respect to a parameter
    
    struct SubstitutionModel * (*clone)( struct SubstitutionModel * );
    
    void (*free)( struct SubstitutionModel * );
	
}SubstitutionModel;

Model * new_SubstitutionModel2( const char* name, SubstitutionModel *sm, Model* simplex  );

void generale_update_freqs( SubstitutionModel *model );

void update_eigen_system( SubstitutionModel *m );

void check_frequencies( const double *freqs, const int dim );

void make_zero_rows( double **q, const int dim );

void normalize_Q( double **m, const double *freqs, const int dim );

Model* new_SubstitutionModel_from_json(json_node* node, Hashtable*hash);

SubstitutionModel * create_substitution_model( const char *name, const modeltype modelname, const datatype dtype, Simplex* freqs );

SubstitutionModel * create_nucleotide_model( const char *name, const modeltype modelname, Simplex* freqs );

SubstitutionModel * create_codon_model( const char *name, const modeltype modelname, unsigned gen_code, Simplex* freqs );

SubstitutionModel * create_aa_model( const char *name, const modeltype modelname, Simplex* freqs );

SubstitutionModel * create_general_model( const char *name, const modeltype modelname, Simplex* freqs );


SubstitutionModel * clone_substitution_model(SubstitutionModel *m);

SubstitutionModel * clone_substitution_model_with(SubstitutionModel *m, const Parameters* rates, Simplex* simplex);

void free_SubstitutionModel( SubstitutionModel *m);

double get_frequency( SubstitutionModel *m, int pos );



#pragma mark -
#pragma mark misc
	
void print_frequencies(FILE *pfile, SubstitutionModel *m);
void print_rates(FILE *pfile, const SubstitutionModel *m);

void bufferize_frequencies(StringBuffer *buffer, SubstitutionModel *m );
void bufferize_codon_frequencies(StringBuffer *buffer, SubstitutionModel *m );
void bufferize_aa_frequencies(StringBuffer *buffer, SubstitutionModel *m );
void bufferize_rates( StringBuffer *buffer, const SubstitutionModel *m );

void compare_model( const SubstitutionModel *m1, const SubstitutionModel *m2 );

SubstitutionModel * SubstitutionModel_factory( const char* model_string, DataType* datatype, Simplex* freqSimplex, const Parameters* rates, const char** assignment );

/*SubstitutionModel * SubstitutionModel_factory( const char* model_string );

SubstitutionModel * SubstitutionModel_nuc_factory_with_values( modeltype modtype, double *rates, double *freqs );

char * SubstitutionModel_get_string( modeltype modtype );*/


void print_P( double *p, int nstate );
void print_P2( double **p, int nstate );
void print_substitution_matrix( SubstitutionModel *sm, double t );

modeltype SubstitutionModel_get_code( const char *model_string );

double tstv_to_kappa( double tstv, const double *freqs );

double kappa_to_tstv( double kappa, const double *freqs );

/*char *SubstitutionModel_stringify( SubstitutionModel *m );
 
 StringBuffer * SubstitutionModel_bufferize( StringBuffer *buffer, SubstitutionModel *m );
 
 void * SubstitutionModel_SML_to_object( ObjectStore *store, SMLNode node );*/

#endif
