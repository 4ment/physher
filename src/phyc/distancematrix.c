/*
 *  distancematrix.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 16/01/2015.
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

#include "distancematrix.h"

#include <ctype.h>
#include <assert.h>
#include <strings.h>

#include "matrix.h"



static bool const AMINOACID_STATES[128] = {
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,	// 0-15
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,	// 16-31
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,	// 32-47
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,	// 48-63
    //	  A    B     C    D    E    F    G    H    I          K    L    M    N
    false,true,false,true,true,true,true,true,true,true,false,true,true,true,true,false,	            // 64-79
    // P Q    R    S    T          V    W    X     Y
    true,true,true,true,true,false,true,true,false,true,false,false,false,false,false,false,	        // 80-95
    //	  a          c    d    e    f    g    h    i          k    l    m    n
    false,true,false,true,true,true,true,true,true,true,false,true,true,true,true,false,	            // 96-111
    // p q    r    s    t          v    w          y
    true,true,true,true,true,false,true,true,false,true,false,false,false,false,false,false		        // 112-127
};

static const int NUCLEOTIDE_STATES[128] = {
    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 0-15
    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 16-31
    //                                          -
    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 32-47
    //                                                ?
    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,	// 48-63
    //	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
    17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 64-79
    //	 p  q  R  S  T  U  V  W  x  Y  z
    16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17,	// 80-95
    //	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
    17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 96-111
    //p q  R  S  T  U  V  W  x  Y  z
    16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17		// 112-127
};

static bool const IS_TRANSITION[4][4] = {
    {false, false, true,  false},
    {false, false, false, true},
    {true,  false, false, false},
    {false, true,  false, false}
};




void _distance_raw(const Sequences *sequences, double **matrix){
	char *seq1 = NULL;
	char *seq2 = NULL;
	int i,j;
	int n1,n2;
	
	int length = 0;
	
	for ( i = 0; i < sequences->size; i++ ) {
		matrix[i][i] = 0.0;
		for ( j = i+1; j < sequences->size; j++ ) {
			seq1 = sequences->seqs[i]->seq;
			seq2 = sequences->seqs[j]->seq;
			double d = 0;
			int n = 0;
			length = sequences->seqs[0]->length;
			
			while(length){
				n1 = NUCLEOTIDE_STATES[*seq1];
				n2 = NUCLEOTIDE_STATES[*seq2];
				if(n1 < 4 && n2 < 4){
					if( n1 != n2 ){
						++d;
					}
					++n;
				}
				++seq1;
				++seq2;
				--length;
			}
			if( n == 0 ){
				matrix[j][i] = matrix[i][j] = 1000.0;
			}
			else {
				matrix[j][i] = matrix[i][j] = d/n;
			}
		}
	}
}

void _distance_raw_patterns(const SitePattern *patterns, double **matrix){
	int i,j;
	int n1,n2;
	
	for ( i = 0; i < patterns->size; i++ ) {
		matrix[i][i] = 0.0;
		for ( j = i+1; j < patterns->size; j++ ) {
			double d = 0;
			int n = 0;
			
			for ( int site = 0; site < patterns->count; site++) {
				n1 = patterns->patterns[site][i];
				n2 = patterns->patterns[site][j];
				if(n1 < 4 && n2 < 4){
					if( n1 != n2 ){
						d += patterns->weights[site];
					}
					n += patterns->weights[site];
				}
			}
			if( n == 0 ){
				matrix[j][i] = matrix[i][j] = 1000.0;
			}
			else {
				matrix[j][i] = matrix[i][j] = d/n;
			}
		}
	}
}

void _distance_raw_float(const Sequences *sequences, float **matrix){
    char *seq1 = NULL;
    char *seq2 = NULL;
   	int i,j;
    int n1,n2;
	
    int length = 0;
	
    for ( i = 0; i < sequences->size; i++ ) {
        matrix[i][i] = 0.f;
        for ( j = i+1; j < sequences->size; j++ ) {
            seq1 = sequences->seqs[i]->seq;
            seq2 = sequences->seqs[j]->seq;
            float d = 0;
            int n = 0;
            length = sequences->seqs[0]->length;
            
            while(length){
                n1 = NUCLEOTIDE_STATES[*seq1];
                n2 = NUCLEOTIDE_STATES[*seq2];
                if(n1 < 4 && n2 < 4){
                    if( n1 != n2 ){
                        ++d;
                    }
                    ++n;
                }
                ++seq1;
                ++seq2;
                --length;
            }
            if( n == 0 ){
                matrix[j][i] = matrix[i][j] = 1000.f;
            }
            else {
                matrix[j][i] = matrix[i][j] = d/n;
            }
        }
    }
}


void _distance_jc69(const Sequences *sequences, double **matrix){
	char *seq1 = NULL;
	char *seq2 = NULL;
	int i,j;
	double d = 0;
	int n = 0;
	double fourThird = 4.f/3.f;
	int length;
	int n1,n2;
	
	for ( i = 0; i < sequences->size; i++ ) {
		matrix[i][i] = 0.0;
		for ( j = i+1; j < sequences->size; j++ ) {
			seq1 = sequences->seqs[i]->seq;
			seq2 = sequences->seqs[j]->seq;
			d = 0.f;
			n = 0;
			
			length = sequences->seqs[0]->length;
			
			while(length){
				n1 = NUCLEOTIDE_STATES[*seq1];
				n2 = NUCLEOTIDE_STATES[*seq2];
				if(n1 < 4 && n2 < 4){
					if( n1 != n2 ){
						++d;
					}
					++n;
				}
				++seq1;
				++seq2;
				--length;
			}
			d /= n;
			if( d >= 0.75 ){
				fprintf(stderr, "Dissimilarity over 0.75 %s %s %f\n", sequences->seqs[i]->name,sequences->seqs[j]->name, d);
				matrix[j][i] = matrix[i][j] = 1000.0;
			}
			else {
				matrix[j][i] = matrix[i][j] = -0.75*log(1.-fourThird*d);
			}
		}
	}
}

void _distance_jc69_patterns(const SitePattern *patterns, double **matrix){
	int i,j;
	double d = 0;
	int n = 0;
	double fourThird = 4.f/3.f;
	int n1,n2;
	
	for ( i = 0; i < patterns->size; i++ ) {
		matrix[i][i] = 0.0;
		for ( j = i+1; j < patterns->size; j++ ) {
			d = 0.f;
			n = 0;
			
			for ( int site = 0; site < patterns->count; site++) {
				n1 = patterns->patterns[site][i];
				n2 = patterns->patterns[site][j];
				if(n1 < 4 && n2 < 4){
					if( n1 != n2 ){
						d += patterns->weights[site];
					}
					n += patterns->weights[site];
				}
			}
			d /= n;
			if( d >= 0.75 ){
				fprintf(stderr, "Dissimilarity over 0.75 %s %s %f\n", patterns->names[i], patterns->names[j], d);
				matrix[j][i] = matrix[i][j] = 1000.0;
			}
			else {
				matrix[j][i] = matrix[i][j] = -0.75*log(1.-fourThird*d);
			}
		}
	}
}

void _distance_jc69_float(const Sequences *sequences, float **matrix){
    char *seq1 = NULL;
    char *seq2 = NULL;
   	int i,j;
    float d = 0;
    int n = 0;
    float fourThird = 4.f/3.f;
    int length;
    int n1,n2;
	
    for ( i = 0; i < sequences->size; i++ ) {
        matrix[i][i] = 0.0;
        for ( j = i+1; j < sequences->size; j++ ) {
            seq1 = sequences->seqs[i]->seq;
            seq2 = sequences->seqs[j]->seq;
            d = 0.f;
            n = 0;
            
            length = sequences->seqs[0]->length;
            
            while(length){
                n1 = NUCLEOTIDE_STATES[*seq1];
                n2 = NUCLEOTIDE_STATES[*seq2];
                if(n1 < 4 && n2 < 4){
                    if( n1 != n2 ){
                        ++d;
                    }
                    ++n;
                }
                ++seq1;
                ++seq2;
                --length;
            }
            d /= n;
            if( d >= 0.75 ){
                fprintf(stderr, "Dissimilarity over 0.75 %s %s %f\n", sequences->seqs[i]->name,sequences->seqs[j]->name, d);
                matrix[j][i] = matrix[i][j] = 1000.0;
            }
            else {
                matrix[j][i] = matrix[i][j] = -0.75f*logf(1.f-fourThird*d);
            }
        }
    }
}

void _distance_k2p(const Sequences *sequences, double **matrix){
    char *seq1 = NULL;
    char *seq2 = NULL;
    double transition, transversion, n;
   	int i,j;
    int length;
    int n1,n2;
    
    for ( i = 0; i < sequences->size; i++ ) {
        matrix[i][i] = 0.0;
        for ( j = i+1; j < sequences->size; j++ ) {
            seq1 = sequences->seqs[i]->seq;
            seq2 = sequences->seqs[j]->seq;
            transition = transversion = n = 0;
            
            length = sequences->seqs[0]->length;
            
            while(length){
                n1 = NUCLEOTIDE_STATES[*seq1];
                n2 = NUCLEOTIDE_STATES[*seq2];
                if(n1 < 4 && n2 < 4){
                    if( IS_TRANSITION[n1][n2]){
                        ++transition;
                    }
                    else if ( n1 != n2 ) {
                        ++transversion;
                    }
                    ++n;
                }
                ++seq1;
                ++seq2;
                --length;
            }
            matrix[j][i] = matrix[i][j] = 0.5 * log(1./(1.-(2.*transition-transversion)/n)) + 0.25*log(1./(1.-2.*transversion/n));
        }
    }
}


void _distance_k2p_patterns(const SitePattern *patterns, double **matrix){
	double transition, transversion, n;
	int i,j;
	int n1,n2;
	
	for ( i = 0; i < patterns->size; i++ ) {
		matrix[i][i] = 0.0;
		for ( j = i+1; j < patterns->size; j++ ) {
			transition = transversion = n = 0;
			
			for ( int site = 0; site < patterns->count; site++) {
				n1 = patterns->patterns[site][i];
				n2 = patterns->patterns[site][j];
				if(n1 < 4 && n2 < 4){
					if( IS_TRANSITION[n1][n2]){
						transition += patterns->weights[site];
					}
					else if ( n1 != n2 ) {
						transversion += patterns->weights[site];
					}
					n += patterns->weights[site];
				}
			}
			matrix[j][i] = matrix[i][j] = 0.5 * log(1./(1.-(2.*transition-transversion)/n)) + 0.25*log(1./(1.-2.*transversion/n));
		}
	}
}

void _distance_k2p_float(const Sequences *sequences, float **matrix){
    char *seq1 = NULL;
    char *seq2 = NULL;
    float transition, transversion, n;
   	int i,j;
    int length;
    int n1,n2;
	
    for ( i = 0; i < sequences->size; i++ ) {
        matrix[i][i] = 0.0;
        for ( j = i+1; j < sequences->size; j++ ) {
            seq1 = sequences->seqs[i]->seq;
            seq2 = sequences->seqs[j]->seq;
            transition = transversion = n = 0;
            
            length = sequences->seqs[0]->length;
            
            while(length){
                n1 = NUCLEOTIDE_STATES[*seq1];
                n2 = NUCLEOTIDE_STATES[*seq2];
                if(n1 < 4 && n2 < 4){
                    if( IS_TRANSITION[n1][n2]){
                        ++transition;
                    }
                    else if ( n1 != n2 ) {
                        ++transversion;
                    }
                    ++n;
                }
                ++seq1;
                ++seq2;
                --length;
            }
            matrix[j][i] = matrix[i][j] = 0.5f * logf(1.f/(1.f-(2.f*transition-transversion)/n)) + 0.25f*logf(1.f/(1.f-2.f*transversion/n));
        }
    }
}

#pragma mark *** Amino acid ***

void _distance_aa_kimura(const Sequences *sequences, double **matrix){
    char *seq1 = NULL;
    char *seq2 = NULL;
    char n1,n2;
   	int i,j;
    double p, d,n;
    
    for ( i = 0; i < sequences->size; i++ ) {
        matrix[i][i] = 0.0;
        for ( j = i+1; j < sequences->size; j++ ) {
            seq1 = sequences->seqs[i]->seq;
            seq2 = sequences->seqs[j]->seq;
            n = d = 0;
            
            for ( int site = 0; site < sequences->seqs[0]->length; site++) {
                n1 = seq1[site];
                n2 = seq2[site];
                
                if ( AMINOACID_STATES[n1] && AMINOACID_STATES[n2] ) {
                    if( n1 != n2 ){
                        d++;
                    }
                    n++;
                }
            }
            p = d/n;
            matrix[j][i] = matrix[i][j] = -log(1.0-p-0.2*p*p);
        }
    }
}

void _distance_aa_kimura_patterns(const SitePattern *patterns, double **matrix){
	int n1,n2;
	int i,j;
	double p, d,n;
	DataType* datatype = patterns->datatype;
	
	for ( i = 0; i < patterns->size; i++ ) {
		matrix[i][i] = 0.0;
		for ( j = i+1; j < patterns->size; j++ ) {
			n = d = 0;
			
			for ( int site = 0; site < patterns->count; site++) {
				n1 = patterns->patterns[site][i];
				n2 = patterns->patterns[site][j];
				if ( n1 < datatype->state_count(datatype) && n2 < datatype->state_count(datatype) ) {
					if( n1 != n2 ){
						d+= patterns->weights[site];
					}
					n += patterns->weights[site];
				}
			}
			p = d/n;
			matrix[j][i] = matrix[i][j] = -log(1.0-p-0.2*p*p);
		}
	}
}

void _distance_aa_kimura_float(const Sequences *sequences, float **matrix){
    char *seq1 = NULL;
    char *seq2 = NULL;
   	int i,j;
    float d,n;
    int length;
	
    for ( i = 0; i < sequences->size; i++ ) {
        matrix[i][i] = 0.0;
        for ( j = i+1; j < sequences->size; j++ ) {
            seq1 = sequences->seqs[i]->seq;
            seq2 = sequences->seqs[j]->seq;
            n = d = 0;
            
            length = sequences->seqs[0]->length;
            
            while(length){
                if ( AMINOACID_STATES[*seq1] && AMINOACID_STATES[*seq2] ) {
                    if( *seq1 != *seq2 ){
                        ++d;
                    }
                    ++n;
                }
                --length;
                ++seq1;
                ++seq2;
            }
            d /= n;
            matrix[j][i] = matrix[i][j] = -logf(1.f-d-0.2f*d*d);
        }
    }
}


void _hamming_distance(const Sequences *sequences, double **matrix){
	char *seq1 = NULL;
	char *seq2 = NULL;
	int i,j;
	
	int length = 0;
	
	for ( i = 0; i < sequences->size; i++ ) {
		matrix[i][i] = 0.f;
		for ( j = i+1; j < sequences->size; j++ ) {
			seq1 = sequences->seqs[i]->seq;
			seq2 = sequences->seqs[j]->seq;
			double d = 0;
			length = sequences->seqs[0]->length;
			
			while(length){
				if( *seq1 != *seq2 ){
					++d;
				}
				
				++seq1;
				++seq2;
				--length;
			}
			matrix[j][i] = matrix[i][j] = d/sequences->seqs[0]->length;
		}
	}
}

void _hamming_distance_patterns(const SitePattern *patterns, double **matrix){
	int i,j;
	int n1,n2;
	
	for ( i = 0; i < patterns->size; i++ ) {
		matrix[i][i] = 0.f;
		for ( j = i+1; j < patterns->size; j++ ) {
			double d = 0;
			
			for ( int site = 0; site < patterns->count; site++) {
				n1 = patterns->patterns[site][i];
				n2 = patterns->patterns[site][j];
				if( n1 != n2 ){
					d += patterns->weights[site];
				}
			}
			matrix[j][i] = matrix[i][j] = d/patterns->count;
		}
	}
}

void _hamming_distance_float(const Sequences *sequences, float **matrix){
    char *seq1 = NULL;
    char *seq2 = NULL;
   	int i,j;
	
    int length = 0;
	
    for ( i = 0; i < sequences->size; i++ ) {
        matrix[i][i] = 0.f;
        for ( j = i+1; j < sequences->size; j++ ) {
            seq1 = sequences->seqs[i]->seq;
            seq2 = sequences->seqs[j]->seq;
            float d = 0;
            length = sequences->seqs[0]->length;
            
            while(length){
                if( *seq1 != *seq2 ){
                    ++d;
                }
                
                ++seq1;
                ++seq2;
                --length;
            }
            matrix[j][i] = matrix[i][j] = d/sequences->seqs[0]->length;
        }
    }
}

double ** Sequences_distance( const Sequences *sequences, distancematrix_model model ){
    
    double **matrix = dmatrix(sequences->size, sequences->size);
    
    if( sequences->datatype->type == DATA_TYPE_AMINO_ACID ){
        _distance_aa_kimura(sequences, matrix);
    }
    else if( sequences->datatype->type == DATA_TYPE_NUCLEOTIDE || sequences->datatype->type == DATA_TYPE_CODON ){
        switch (model) {
            case DISTANCE_MATRIX_UNCORRECTED:
                _distance_raw(sequences, matrix);
                break;
            case DISTANCE_MATRIX_JC69:
                _distance_jc69(sequences, matrix);
                break;
            case DISTANCE_MATRIX_K2P:
                _distance_k2p(sequences, matrix);
                break;
            default:
                assert(0);
        }
    }
    else {
        _hamming_distance(sequences, matrix);
    }
    return matrix;
}

float ** Sequences_distance_float( const Sequences *sequences, distancematrix_model model ){
    
    float **matrix = fmatrix(sequences->size, sequences->size);
    
    if( sequences->datatype->type == DATA_TYPE_AMINO_ACID ){
        _distance_aa_kimura_float(sequences, matrix);
    }
    else if( sequences->datatype->type == DATA_TYPE_NUCLEOTIDE || sequences->datatype->type == DATA_TYPE_CODON ){
        
        switch (model) {
            case DISTANCE_MATRIX_UNCORRECTED:
                _distance_raw_float(sequences, matrix);
                break;
            case DISTANCE_MATRIX_JC69:
                _distance_jc69_float(sequences, matrix);
                break;
            case DISTANCE_MATRIX_K2P:
                _distance_k2p_float(sequences, matrix);
                break;
            default:
                assert(0);
        }
    }
    else {
        _hamming_distance_float(sequences, matrix);
    }
    return matrix;
}

double ** SitePattern_distance( const SitePattern *patterns, distancematrix_model model ){
	
	double **matrix = dmatrix(patterns->size, patterns->size);
	
	if( patterns->datatype->type == DATA_TYPE_AMINO_ACID ){
		_distance_aa_kimura_patterns(patterns, matrix);
	}
	else if( patterns->datatype->type == DATA_TYPE_NUCLEOTIDE || patterns->datatype->type == DATA_TYPE_CODON ){
		switch (model) {
			case DISTANCE_MATRIX_UNCORRECTED:
				_distance_raw_patterns(patterns, matrix);
				break;
			case DISTANCE_MATRIX_JC69:
				_distance_jc69_patterns(patterns, matrix);
				break;
			case DISTANCE_MATRIX_K2P:
				_distance_k2p_patterns(patterns, matrix);
				break;
			default:
				assert(0);
		}
	}
	else {
		_hamming_distance_patterns(patterns, matrix);
	}
	return matrix;
}

Matrix* create_DistanceMatrix_from_json( json_node* node, Hashtable* hash ){
	SitePattern* patterns = NULL;
	json_node* patterns_node = get_json_node(node, "sitepattern");
	char* model = get_json_node_value_string(node, "model");
	if(patterns_node->node_type == MJSON_STRING){
		char* ref = (char*)patterns_node->value;
		// check it starts with a &
		patterns = Hashtable_get(hash, ref+1);
		patterns->ref_count++;
	}
	else{
		char* id = get_json_node_value_string(patterns_node, "id");
		patterns = new_SitePattern_from_json(patterns_node, hash);
		Hashtable_add(hash, id, patterns);
	}
	distancematrix_model modelt = DISTANCE_MATRIX_UNCORRECTED;
	if (strcasecmp(model, "uncorrected")) {
		modelt = DISTANCE_MATRIX_UNCORRECTED;
	}
	else if (strcasecmp(model, "jc69")) {
		modelt = DISTANCE_MATRIX_JC69;
	}
	else if (strcasecmp(model, "k2p")) {
		modelt = DISTANCE_MATRIX_K2P;
	}
	else if (strcasecmp(model, "kimura")) {
		modelt = DISTANCE_MATRIX_KIMURA;
	}
	double** matrix = SitePattern_distance(patterns, modelt);
	Matrix * mat = new_Matrix_from(matrix, patterns->size, patterns->size);
	mat->rowNames = malloc(patterns->size*sizeof(char*));
	for (int i = 0; i < patterns->size; i++) {
		mat->rowNames[i] = String_clone(patterns->names[i]);
	}
	free_SitePattern(patterns);
	return mat;
}
