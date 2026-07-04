// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include "utils.h"

#include <string.h>
#include <stdint.h>

#include "datatype.h"
#include "parameters.h"

//typedef enum distance_t{DISTANCE_RAW, DISTANCE_JC69}distance_t;

typedef struct Sequence{
    char *name;
	char *seq;
    int length;
} Sequence;

typedef struct Sequences{
	Sequence **seqs;
    int size;
	int capacity;
	int length;
	bool aligned;
    
	DataType *datatype;
} Sequences;


#pragma mark -
#pragma mark Sequences

Sequences * new_Sequences( unsigned capacity );

void free_Sequences( Sequences *a );

Sequences * Sequences_clone( const Sequences * sequences );

void Sequences_add(Sequences *aln, Sequence *sequence );

void Sequences_delete( Sequences *sequences, int index );

void Sequences_pack( Sequences *seqs );

void Sequences_pack2( Sequences *sequences );

bool Sequences_concatenate( Sequences *seqs1, const Sequences *seqs2 );

int Sequences_get_index( Sequences *seqs, const char *name );

void Sequences_sort_from_ivector( Sequences *seqs, int *s, int size );

Sequences ** Sequences_split3(Sequences *sequences );

Sequences ** Sequences_split001(Sequences *sequences );

#pragma mark -
#pragma mark Sequence


Sequence * new_Sequence( const char *name, const char *seq );

Sequence * new_Sequence2( const char *name, int capacity );

void free_Sequence( Sequence *s );

Sequence * Sequence_clone( const Sequence *sequence );



void empirical_frequencies( const Sequences *sequences, double *freqs );

void empirical_generic_frequencies( const Sequences *sequences, double *freqs );

void empirical_nuc_frequencies( const Sequences *sequences, double *freqs );

void calculate_nuc_freq_codon( const Sequences *sequences, double *nuc_freq );

void calculate_F1x4( const Sequences *sequences, double *codon_freq );

void calculate_F3x4( const Sequences *sequences, double *codon_freq );

void calculate_F1x4_from_nuc_freq( unsigned gen_code, const double *nuc_freq, double *codon_freq );

void calculate_F3x4_from_nuc_freq( unsigned gen_code, const double *nuc_freq, double *codon_freq );

double Sequence_kappa_empirical( const Sequences *sequences, const double *freqs );

double Sequence_tstv_empirical( const Sequences *sequences );

void Sequence_rel_rates_empirical( const Sequences *sequences, double *rel );

Model* new_Alignment_from_json(json_node* node, Hashtable* hash);
	
#endif
