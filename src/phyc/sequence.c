/*
 *  sequence.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/8/10.
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


#include "sequence.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "matrix.h"
#include "utils.h"
#include "mstring.h"
#include "geneticcode.h"
#include "sitepattern.h"
#include "random.h"


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
    //	 p  q  R  S  T  U  V  W  x  Y  z
    16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17		// 112-127
};

static bool const IS_TRANSITION[4][4] = {
    {false, false, true,  false},
    {false, false, false, true},
    {true,  false, false, false},
    {false, true,  false, false}
};

#pragma mark -
#pragma mark Sequences

Sequences * new_Sequences( unsigned capacity ){
	Sequences *sequences = (Sequences*)malloc(sizeof(Sequences));
    assert(sequences);
	sequences->seqs = (Sequence**)calloc(capacity, sizeof(Sequence*));
    assert(sequences->seqs);
	sequences->aligned = true;
    sequences->datatype = NULL;//nucleotide_datatype();
	sequences->length = 0;
	sequences->size = 0;
	sequences->capacity = capacity;
	return sequences;
}



void free_Sequences( Sequences *a ){
	for (int i = 0; i < a->size; i++) {
		free_Sequence( a->seqs[i]);
	}
	free(a->seqs);
    if(a->datatype != NULL)free_DataType(a->datatype);
	free(a);
	
}

Sequences * Sequences_clone( const Sequences * sequences ){
    assert(sequences->datatype);
	Sequences *seqs = new_Sequences(sequences->size);
	for ( int i = 0; i < sequences->size; i++ ) {
		Sequences_add(seqs, sequences->seqs[i]);
	}
	seqs->datatype = clone_DataType(sequences->datatype);
	seqs->aligned = sequences->aligned;
	return seqs;
}

void Sequences_add( Sequences *sequences, Sequence *sequence ){
	if( sequences->size > sequences->capacity-1 ){
		sequences->capacity *= 2;
		sequences->seqs = realloc( sequences->seqs, sequences->capacity * sizeof(Sequence *));
	}
	sequences->seqs[sequences->size] = sequence;
	
	if ( sequences->aligned ) {
		if( sequences->size == 0 ){
			sequences->length = sequence->length;
		}
	}
	else {
		sequences->length = sequence->length;
	}
	
	sequences->size++;
	
}

// free a sequence and move all the sequences after it
void Sequences_delete( Sequences *sequences, int index ){
	free_Sequence(sequences->seqs[index]);
	sequences->seqs[index] = NULL;
    while ( index < sequences->size-1 ) {
        sequences->seqs[index] = sequences->seqs[index+1];
        index++;
    }
    sequences->seqs[index] = NULL;
	sequences->size--;
}

void Sequences_pack2( Sequences *sequences ){
    
    while ( sequences->seqs[sequences->size-1] == NULL ) {
        sequences->size--;
    }
    
    for( int i = 0; i < sequences->size; i++ ){
        if( sequences->seqs[i] == NULL ){
            int j = i;
            while ( sequences->seqs[j] == NULL && j != sequences->size ) {
                j++;
            }
            if( j != sequences->size ){
                sequences->seqs[i] = sequences->seqs[j];
                sequences->seqs[j] = NULL;
            }
        }
    }
    while ( sequences->seqs[sequences->size-1] == NULL ) {
        sequences->size--;
    }
//    int i = 0;
//    while ( i != sequences->size-1 ) {
//        if( sequences->seqs[i] == NULL ){
//            int j = i;
////            while ( j != sequences->size-2 ) {
////                sequences->seqs[j] = sequences->seqs[j+1];
////                j++;
////            }
//            
//            while ( sequences->seqs[j] == NULL ) {
//                j++;
//            }
//            sequences->seqs[i] = sequences->seqs[j];
//            sequences->size--;
//            sequences->seqs[sequences->size] = NULL;
//        }
//        else {
//            i++;
//        }
//    }
    //printf("size %d\n",sequences->size);
}

void Sequences_pack( Sequences *seqs ){
	if ( seqs->size != seqs->capacity ) seqs->seqs = realloc( seqs->seqs, seqs->size * sizeof(Sequence *));
	seqs->capacity = seqs->size;
}


bool Sequences_concatenate( Sequences *seqs1, const Sequences *seqs2 ){
	if ( seqs1->size != seqs2->size ) {
		return false;
	}
	int *indexes = ivector(seqs1->size);
	int j = 0;
	for ( int i = 0; i < seqs1->size; i++ ) {
		for ( j = 0; j < seqs1->size; j++ ) {
			if( strcmp( seqs1->seqs[i]->name, seqs2->seqs[j]->name) == 0 ){
				indexes[i] = j;
				break;
			}
		}
		if ( j == seqs1->size ) {
			free(indexes);
			return false;
		}
	}
	
	for ( int i = 0; i < seqs1->size; i++ ) {
		String_append_string( seqs1->seqs[i]->seq, seqs2->seqs[ indexes[i] ]->seq );
		seqs1->length += seqs2->length;
	}
	seqs1->size *= 2;
	free(indexes);
	return true;
}

int Sequences_get_index( Sequences *seqs, const char *name ){
    for ( int i = 0; i < seqs->size; i++ ) {
        if( strcmp( seqs->seqs[i]->name, name) == 0 ){
            return i;
        }
    }
    return -1;
}

void Sequences_sort_from_ivector( Sequences *seqs, int *s, int size ){
	bool done = false;
    Sequence *temp = NULL;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( s[i] > s[i+1] ) {
				done = false;
				swap_int( &s[i], &s[i+1] );
                temp = seqs->seqs[i];
                seqs->seqs[i] = seqs->seqs[i+1];
                seqs->seqs[i+1] = temp;
			}
		}
		size--;
	}
}

#pragma mark -
#pragma mark Sequence

Sequence * new_Sequence( const char *name, const char *seq ){
	Sequence *sequence = (Sequence*)malloc(sizeof(Sequence));
	assert(sequence);
	sequence->name   = String_clone(name);
	sequence->seq    = String_clone(seq);
	sequence->length = strlen(seq);
	return sequence;
}

// parameter length is the length of the sequence and does not include the final char '\0'
Sequence * new_Sequence2( const char *name, int capacity ){
	Sequence *sequence = (Sequence*)malloc(sizeof(Sequence));
	assert(sequence);
	sequence->name   = String_clone(name);
	sequence->seq    = String_new(capacity+1);
	sequence->length = 0;
	for ( int i = 0; i <= capacity; i++ ) {
		sequence->seq[i] = '\0';
	}
	return sequence;
}

void free_Sequence( Sequence *s ){
	free(s->name);
	free(s->seq);
	free(s);
}

Sequence * Sequence_clone( const Sequence *sequence ){
	return new_Sequence(sequence->name, sequence->seq);
}

void calculate_F3x4( const Sequences *sequences, double *codon_freq ){
	double *nuc_freq = dvector(12);
	calculate_nuc_freq_codon( sequences, nuc_freq);
	calculate_F3x4_from_nuc_freq(sequences->datatype->genetic_code, nuc_freq, codon_freq);
	free(nuc_freq);
}

void calculate_F1x4( const Sequences *sequences, double *codon_freq ){
	double *nuc_freq = dvector(4);
	empirical_nuc_frequencies(sequences, nuc_freq);
	calculate_F1x4_from_nuc_freq(sequences->datatype->genetic_code, nuc_freq, codon_freq);
	free(nuc_freq);
}

void calculate_nuc_freq_codon( const Sequences *sequences, double *nuc_freq ){
	double *pos_count = dvector(3);
	int i = 0;
	Sequence *seq = NULL;
    DataType *nuctype = nucleotide_datatype();//SINGLETON_DATATYPE_NUCLEOTIDE;
	
	for ( i = 0; i < sequences->size; i++ ) {
		seq = sequences->seqs[i];
		for ( int site = 0; site < seq->length; site++) {
			int j = nuctype->encoding(nuctype, seq->seq[site]);
            
			if ( j < 4 ){
                int p = site % 3;
                pos_count[p]++;
                nuc_freq[p*4+j]++;
            }
		}
	}
	for ( i = 0; i < 12; i++ ) {
		assert(nuc_freq[i] > 0);
		nuc_freq[i] /= pos_count[i/4];
	}
	free(pos_count);
	
	//	fprintf(stderr, "\n\n");
	//	for ( i = 0; i < 12; i+=4 ) {
	//		fprintf(stderr, "pos %d %f %f %f %f\n",i/4+1,nuc_freq[i],nuc_freq[i+1],nuc_freq[i+2],nuc_freq[i+3]);
	//	}
}

void empirical_frequencies( const Sequences *sequences, double *freqs ){
    DataType *dataType = sequences->datatype;
    assert(dataType);
    int count = 0;
    int i = 0;
    int nstate  = dataType->state_count(dataType);
    
    Sequence *seq = NULL;
    memset(freqs, 0, nstate*sizeof(double));
    
    for ( i = 0; i < sequences->size; i++ ) {
        seq = sequences->seqs[i];
        for ( int site = 0; site < seq->length; site++) {
            int encoding = dataType->encoding(dataType, seq->seq[site]);
            
            if ( encoding < nstate ){
                count++;
                freqs[encoding]++;
            }
        }
    }
    for ( i = 0; i < nstate; i++ ) {
        freqs[i] /= count;
    }
}

void empirical_generic_frequencies( const Sequences *sequences, double *freqs ){
    DataType *dataType = sequences->datatype;
    assert(dataType);
    int count = 0;
    int i = 0;
    int nstate  = dataType->state_count(dataType);
    char *state = cvector(sequences->length+1);
    state[sequences->length] = '\0';
    
    Sequence *seq = NULL;
    memset(freqs, 0, nstate*sizeof(double));
    
    for ( i = 0; i < sequences->size; i++ ) {
        seq = sequences->seqs[i];
        for ( int site = 0; site < seq->length; site++) {
            memcpy(state, seq->seq, seq->length*sizeof(char));
            int encoding = dataType->encoding_string(dataType, state);
            
            if ( encoding < nstate ){
                count++;
                freqs[encoding]++;
            }
        }
    }
    for ( i = 0; i < nstate; i++ ) {
        assert(freqs[i] > 0.0);
        freqs[i] /= count;
    }
    free(state);
}

void empirical_nuc_frequencies( const Sequences *sequences, double *freqs ){
    int count = 0;
    int i = 0;
    int nstate  = sequences->datatype->state_count(sequences->datatype);
    DataType *nuctype = nucleotide_datatype();//SINGLETON_DATATYPE_NUCLEOTIDE;
    
    Sequence *seq = NULL;
    memset(freqs, 0, nstate*sizeof(double));
    
    for ( i = 0; i < sequences->size; i++ ) {
        seq = sequences->seqs[i];
        for ( int site = 0; site < seq->length; site++) {
            int encoding = nuctype->encoding(nuctype, seq->seq[site]);
            
            if ( encoding < nstate ){
                count++;
                freqs[encoding]++;
            }
        }
    }
    for ( i = 0; i < nstate; i++ ) {
        assert(freqs[i] > 0.0);
        freqs[i] /= count;
    }
}


void calculate_F3x4_from_nuc_freq( unsigned gen_code, const double *nuc_freq, double *codon_freq ){	
	int stops = 0;
	double total = 0.0;
	int n1,n2,n3;
	for( n1 = 0; n1 < 4; n1++ ){
		for( n2 = 0; n2 < 4; n2++ ){
			for( n3 = 0; n3 < 4; n3++ ){
				if( GENETIC_CODE_TABLES[gen_code][n1*16+n2*4+n3] != '*' ) {
					codon_freq[n1*16+n2*4+n3 - stops] = nuc_freq[n1] * nuc_freq[4+n2] * nuc_freq[8+n3];
					total += codon_freq[n1*16+n2*4+n3 - stops];
				}
				else stops++;
			}
		}
	}
	
	for( int i = 0; i < (64-stops); i++ ){
		codon_freq[i] /= total;
		//fprintf(stdout, "%d f %f\n", i, codon_freq[i] );
	}
	
//	int ii = 0;
//	for( int i = 0; i < 64; i++ ){
//		if( GENETIC_CODE_TABLES[gen_code][i] == '*' ) continue;
//		fprintf(stdout, "%d %s %c f %f\n", ii, CODON_TRIPLETS[i], GENETIC_CODE_TABLES[gen_code][i], codon_freq[ii] );
//		ii++;
//	}
}


void calculate_F1x4_from_nuc_freq( unsigned gen_code, const double *nuc_freq, double *codon_freq ){	
	int stops = 0;
	double total = 0.0;
	int n1,n2,n3;
	for( n1 = 0; n1 < 4; n1++ ){
		for( n2 = 0; n2 < 4; n2++ ){
			for( n3 = 0; n3 < 4; n3++ ){
				if( GENETIC_CODE_TABLES[gen_code][n1*16+n2*4+n3] != '*' ) {
					codon_freq[n1*16+n2*4+n3 - stops] = nuc_freq[n1] * nuc_freq[n2] * nuc_freq[n3];
					total += codon_freq[n1*16+n2*4+n3 - stops];
				}
				else stops++;
			}
		}
	}
	
	for( int i = 0; i < (64-stops); i++ ){
		codon_freq[i] /= total;
		//fprintf(stdout, "%d f %f\n", i, codon_freq[i] );
	}	
}


static double _tstv_to_kappa( double tstv, const double *freqs ){
    return (tstv*(freqs[0]+freqs[2])*(freqs[1]+freqs[3]))/((freqs[0]*freqs[2])+(freqs[1]*freqs[3]));
}

// freqs: A C G T
double Sequence_kappa_empirical( const Sequences *sequences, const double *freqs ){
    double tstv = Sequence_tstv_empirical(sequences);
	return _tstv_to_kappa(tstv,freqs);
}


double Sequence_tstv_empirical( const Sequences *sequences ){
	char *seq1 = NULL;
	char *seq2 = NULL;
	int transition = 0;
	int transversion = 0;
	int i,j;
    int n1,n2;
    int length;
    
	for ( i = 0; i < sequences->size; i++ ) {
		for ( j = i+1; j < sequences->size; j++ ) {
			seq1 = sequences->seqs[i]->seq;
			seq2 = sequences->seqs[j]->seq;
            
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
                }
                ++seq1;
                ++seq2;
                --length;
			}
		}
	}
	return (double)transition/transversion;
}

//    A C G T
// A  * a b c
// C  a * d e
// G  b d * 1
// T  c e 1 *

// A C G T
// 0 1 2 3

static int const INDEXES[4][4] = {
    {6, 0, 1, 2},
    {0, 6, 3, 4},
    {1, 3, 6, 5},
    {2, 4, 5, 6}
};

// rel: relative rates of length 5
void Sequence_rel_rates_empirical( const Sequences *sequences, double *rel ){
    char *seq1 = NULL;
    char *seq2 = NULL;
    int n1,n2;
    int i,j,site;
    double rel2[6];
    
    memset(rel2, 0, sizeof(double)*6);
    for ( i = 0; i < sequences->size; i++ ) {
        for ( j = i+1; j < sequences->size; j++ ) {
            seq1 = sequences->seqs[i]->seq;
            seq2 = sequences->seqs[j]->seq;
            for ( site = 0; site < sequences->seqs[0]->length; site++) {
                n1 = NUCLEOTIDE_STATES[seq1[site]];
                n2 = NUCLEOTIDE_STATES[seq2[site]];
                
                if( n1 != n2 && n1 < 4 && n2 < 4 ){
                    rel2[INDEXES[n1][n2]]++;
                }
            }
        }
    }
    
    for ( int i = 0; i < 5; i++ ) {
        rel[i] = rel2[i]/rel2[5];
    }
}

void Sequence_rel_rates_empirical2( const Sequences *sequences, double *rel ){
    char *seq1 = NULL;
    char *seq2 = NULL;
    char n1,n2;
    double gt = 0;
    int i,j,site;
    
    //double *rel = dvector(6);
    memset(rel, 0, sizeof(double)*5);
    for ( i = 0; i < sequences->size; i++ ) {
        for ( j = i+1; j < sequences->size; j++ ) {
            seq1 = sequences->seqs[i]->seq;
            seq2 = sequences->seqs[j]->seq;
            for ( site = 0; site < sequences->seqs[0]->length; site++) {
                n1 = (seq1[site] == 'U' ? 'T' : seq1[site]);
                n2 = (seq2[site] == 'U' ? 'T' : seq2[site]);
                
                if( n1 == '-' || n2 == '-' ){
                    
                }
                else if (   (n1 == 'A' && n2 == 'C') || (n1 == 'C' && n2 == 'A') ){     // A C
                    rel[0]++;
                }
                else if (   (n1 == 'A' && n2 == 'G') || (n1 == 'G' && n2 == 'A') ){     // A G
                    rel[1]++;
                }
                else if (   (n1 == 'A' && n2 == 'T') || (n1 == 'T' && n2 == 'A') ){     // A T
                    rel[2]++;
                }
                else if (   (n1 == 'C' && n2 == 'G') || (n1 == 'G' && n2 == 'C') ){     // C G
                    rel[3]++;
                }
                else if (   (n1 == 'C' && n2 == 'T') || (n1 == 'T' && n2 == 'C') ){     // C T
                    rel[4]++;
                }
                else if (   (n1 == 'G' && n2 == 'T') || (n1 == 'T' && n2 == 'G') ){     // G T
                    gt++;
                }
                
            }
        }
    }
    
    //	fprintf(stderr, "Rel rates");
    //	for ( int i = 0; i < 5; i++ ) {
    //		fprintf(stderr, " %f", rel[i]/rel[5]);
    //	}
    //	fprintf(stderr, "\n");
    //    free(rel);
    for ( int i = 0; i < 5; i++ ) {
        rel[i] /= gt;
    }
}


bool Sequence_isAligned( const Sequences *seqs ){
	size_t l = seqs->seqs[0]->length;
	
	for ( int i = 1; i < seqs->size; i++ ) {
		if ( l != seqs->seqs[i]->length ) {
			//fprintf(stderr, "Sequences have different length\n");
			//exit(1);
			return false;
		}
	}
	return true;
}

void Sequence_remove_empty_columns( Sequences *seqs ){
	for ( int i = seqs->length-1; i >= 0; i-- ) {
		int count = 0;
		for ( int j = 0; j < seqs->size; j++ ) {
			if ( seqs->seqs[j]->seq[i] == '-' || seqs->seqs[j]->seq[i] == ' ' || seqs->seqs[j]->seq[i] == '?' ) {
				count++;
			}
		}
		if ( count == seqs->size ) {
			// Empty column at the end of alignment
			if ( i == seqs->length-1 ) {
				for ( int j = 0; j < seqs->size; j++ ) {
					seqs->seqs[j]->seq[i] = '\0';
					seqs->seqs[j]->length--;
				}
			}
			else {
				for ( int j = 0; j < seqs->size; j++ ) {
					memmove(&seqs->seqs[j]->seq[i], &seqs->seqs[j]->seq[i+1], (seqs->length-i)*sizeof(char)  ); // also move \0
					seqs->seqs[j]->length--;
				}
			}
			seqs->length--;
		}
	}
}


Sequences ** Sequences_split3(Sequences *sequences ){
    DataType *dataType = sequences->datatype;
    assert(dataType);
    Sequences ** seqs = (Sequences**)malloc(sizeof(Sequences*)*3);
    assert(seqs);
    for ( int i = 0; i < 3; i++ ) {
        seqs[i] = new_Sequences(sequences->size);
        seqs[i]->datatype = clone_DataType(dataType);
        
    }
    StringBuffer *buffer1 = new_StringBuffer(sequences->length/3);
    StringBuffer *buffer2 = new_StringBuffer(sequences->length/3);
    StringBuffer *buffer3 = new_StringBuffer(sequences->length/3);
    
    for ( int i = 0; i < sequences->size; i++ ) {
        for ( int j = 0; j < sequences->length; j+=3 ) {
            StringBuffer_append_char(buffer1, sequences->seqs[i]->seq[j] );
            StringBuffer_append_char(buffer2, sequences->seqs[i]->seq[j+1] );
            StringBuffer_append_char(buffer3, sequences->seqs[i]->seq[j+2] );
        }
        
        Sequences_add(seqs[0], new_Sequence(sequences->seqs[i]->name, buffer1->c));
        Sequences_add(seqs[1], new_Sequence(sequences->seqs[i]->name, buffer2->c));
        Sequences_add(seqs[2], new_Sequence(sequences->seqs[i]->name, buffer3->c));
        
        StringBuffer_empty(buffer1);
        StringBuffer_empty(buffer2);
        StringBuffer_empty(buffer3);
    }
    
    free_StringBuffer(buffer1);
    free_StringBuffer(buffer2);
    free_StringBuffer(buffer3);
    
    return seqs;
    
}

// split sequences codon position 1 and 2 from position 3
Sequences ** Sequences_split001(Sequences *sequences ){
    DataType *dataType = sequences->datatype;
    assert(dataType);
    Sequences ** seqs = (Sequences**)malloc(sizeof(Sequences*)*2);
    assert(seqs);
    for ( int i = 0; i < 2; i++ ) {
        seqs[i] = new_Sequences(sequences->size);
        seqs[i]->datatype = clone_DataType(dataType);
        
    }
    StringBuffer *buffer1 = new_StringBuffer(sequences->length);
    StringBuffer *buffer2 = new_StringBuffer(sequences->length/3);
    
    for ( int i = 0; i < sequences->size; i++ ) {
        for ( int j = 0; j < sequences->length; j+=3 ) {
            StringBuffer_append_char(buffer1, sequences->seqs[i]->seq[j] );
            StringBuffer_append_char(buffer1, sequences->seqs[i]->seq[j+1] );
            StringBuffer_append_char(buffer2, sequences->seqs[i]->seq[j+2] );
        }
        
        Sequences_add(seqs[0], new_Sequence(sequences->seqs[i]->name, buffer1->c));
        Sequences_add(seqs[1], new_Sequence(sequences->seqs[i]->name, buffer2->c));
        
        StringBuffer_empty(buffer1);
        StringBuffer_empty(buffer2);
    }
    
    free_StringBuffer(buffer1);
    free_StringBuffer(buffer2);
    
    return seqs;
    
}

