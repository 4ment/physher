/*
 *  datatype.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 19/09/2014.
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

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "datatype.h"
#include "geneticcode.h"
#include "mstring.h"

static int  _state_count( DataType *datatype);


static int  _nucleotide_encoding( DataType *datatype, char nuc);
static char _nucleotide_state( DataType *datatype, int encoding);
static int _nucleotide_encoding_string( DataType *datatype, const char *nuc);
static const char * _nucleotide_state_string( const DataType *datatype, int encoding);

static int  _aa_encoding( DataType *datatype, char nuc);
static char _aa_state( DataType *datatype, int encoding);
static int _aa_encoding_string( DataType *datatype, const char *aa);
static const char * _aa_state_string( const DataType *datatype, int encoding);

static int  _codon_encoding( DataType *datatype, char nuc);
static char _codon_state( DataType *datatype, int encoding);
static int _codon_encoding_string( DataType *datatype, const char *codon);
static const char * _codon_state_string( const DataType *datatype, int encoding);

static DataType SINGLETON_DATATYPE_NUCLEOTIDE = {"Nucleotide", DATA_TYPE_NUCLEOTIDE,  4, 1, 0, _nucleotide_encoding, _nucleotide_state, _nucleotide_encoding_string, _nucleotide_state_string, _state_count, 0};
static DataType SINGLETON_DATATYPE_AMINO_ACID = {"Amino Acid", DATA_TYPE_AMINO_ACID, 20, 1, 0, _aa_encoding,         _aa_state,         _aa_encoding_string,         _aa_state_string        , _state_count, 0};
static DataType SINGLETON_DATATYPE_CODON      = {"Codon",      DATA_TYPE_CODON,      60, 3, 0, _codon_encoding,      _codon_state,      _codon_encoding_string,      _codon_state_string     , _state_count,  0};

static char AMINO_ACIDS[26] = "ACDEFGHIKLMNPQRSTVWYBZX*?-";
static const char *AMINO_ACIDS_STRING[26] = {"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","B","Z","X","*","?","-"};
static int AMINO_ACID_STATES[128] = {
    25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,	// 0-15
    25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,	// 16-31
    //                                 *        -
    25,25,25,25,25,25,25,25,25,25,23,25,25,25,25,25,	// 32-47
    //                                                ?
    25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,24,	// 48-63
    //		A  B  C  D  E  F  G  H  I  j  K  L  M  N  o
    25, 0,20, 1, 2, 3, 4, 5, 6, 7,24, 8, 9,10,11,24,	// 64-79
    //	 P  Q  R  S  T  u  V  W  X  Y  Z
    12,13,14,15,16,24,17,18,22,19,21,25,25,25,25,25,	// 80-95
    //		A  B  C  D  E  F  G  H  I  j  K  L  M  N  o
    25, 0,20, 1, 2, 3, 4, 5, 6, 7,24, 8, 9,10,11,24,	// 96-111
    //	 P  Q  R  S  T  u  V  W  X  Y  Z
    12,13,14,15,16,24,17,18,22,19,21,25,25,25,25,25		// 112-127
};

static const char NUCLEOTIDES[18] = "ACGTUKMRSWYBDHVN?-";
static const char *NUCLEOTIDES_STRING[18] = {"A","C","G","T","U","K","M","R","S","W","Y","B","D","H","V","N","?","-"};
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

int  _state_count( DataType *datatype){
    return datatype->stateCount;
}

const char * _state_string( const DataType *datatype, int encoding){
    if( encoding < datatype->stateCount ){
        return datatype->states[encoding];
    }
    return "?";
}

DataType* new_DataType_from_json(json_node* node, Hashtable* hash){
	DataType* datatype = NULL;
	if(node->node_type == MJSON_OBJECT){
		json_node* states_node = get_json_node(node, "states");
		char* name = get_json_node_value_string(node, "id");
		size_t count = states_node->child_count;
		const char** states = malloc(sizeof(char*)*count);
		for (size_t i = 0; i < count; i++) {
			json_node* child = states_node->children[i];
			states[i] = String_clone((char*)child->value);
		}
		datatype = new_GenericDataType(name, count, states);
		for (size_t i = 0; i < count; i++) {
			free(states[i]);
		}
		free(states);
	}
	else if(node->node_type == MJSON_STRING){
		char* ref = (char*)node->value;
		if (ref[0] == '&') {
			datatype = Hashtable_get(hash, ref+1);
			datatype = clone_DataType(datatype);
		}
		else if(strcasecmp("nucleotide", ref) == 0){
			datatype = new_NucleotideDataType();
		}
		else if(strcasecmp("codon", ref) == 0){
			json_node* code_node = get_json_node(node, "code");
			datatype = new_CodonDataType(atoi((char*)code_node->value));
		}
		else if(strcasecmp("aa", ref) == 0){
			datatype = new_AminoAcidDataType();
		}
		else{
			printf("Datatype not recognixzed %s\n", ref);
			exit(1);
		}
	}
	datatype->ref_count = 1;
	return datatype;
}

DataType * clone_DataType(const DataType *dataType){
	printf("clone datatype\n");
    DataType *newDataType = malloc(sizeof(DataType));
    assert(dataType);
    newDataType->name = String_clone(dataType->name);
    newDataType->type = dataType->type;
    newDataType->stateCount = dataType->stateCount;
    newDataType->symbolLength = dataType->symbolLength;
    newDataType->genetic_code = dataType->genetic_code;
    newDataType->states = malloc(sizeof(char*)*dataType->stateCount);
    assert(dataType->states);
    for ( int i = 0; i < dataType->stateCount; i++ ) {
        newDataType->states[i] = String_clone(dataType->states[i]);
    }
    newDataType->encoding = dataType->encoding;
    newDataType->state    = dataType->state;
    
    newDataType->encoding_string = dataType->encoding_string;
    newDataType->state_string    = dataType->state_string;
    newDataType->state_count = dataType->state_count;
	newDataType->ref_count = 1;
    return newDataType;
}

#pragma mark-
#pragma mark Generic

static int _encoding_string( DataType *datatype, const char *nuc){
    int i = 0;
    for (; i < datatype->stateCount; i++ ) {
        if( strcmp(datatype->states[i], nuc) == 0){
            return i;
        }
    }
    return i;
}

static int _encoding( DataType *datatype, char nuc){
    int i = 0;
    for (; i < datatype->stateCount; i++ ) {
        if( datatype->states[i][0] == nuc ){
            return i;
        }
    }
    exit(1);
    return i;
}

DataType *new_GenericDataType( const char* name, size_t count, const char **states){
	printf("new datatype\n");
    DataType *dataType = malloc(sizeof(DataType));
    assert(dataType);
    dataType->name = String_clone(name);
    dataType->type = DATA_TYPE_GENERIC;
    dataType->genetic_code = 0;
    dataType->stateCount = count;
    dataType->symbolLength = strlen(states[0]);
    dataType->states = malloc(count*sizeof(char*));
    assert(dataType->states);
    for ( int i = 0; i < count; i++ ) {
        dataType->states[i] = String_clone(states[i]);
    }
    dataType->encoding = NULL;
    dataType->state = NULL;
    
    dataType->encoding_string = _encoding_string;
    dataType->state_string    = _state_string;
    
    dataType->state_count = _state_count;
    
    return dataType;
}


void free_DataType( DataType *dataType ){
	if(dataType->ref_count == 1){
//		printf("free datatype\n");
		if(dataType->states != NULL ){
			for ( int i = 0; i < dataType->stateCount; i++ ) {
				free(dataType->states[i]);
			}
			free(dataType->states);
		}
		free(dataType->name);
		free(dataType);
	}
	else{
		dataType->ref_count--;
	}
}



#pragma mark-
#pragma mark Nucleotide

DataType *new_NucleotideDataType(){
    DataType *dataType = (DataType*)malloc(sizeof(DataType));
    assert(dataType);
    dataType->name = String_clone("Nucleotide");
    dataType->type = DATA_TYPE_NUCLEOTIDE;
    dataType->genetic_code = 0;
    dataType->stateCount = 4;
    dataType->symbolLength = 1;
    dataType->states = malloc(4*sizeof(char*));
    assert(dataType->states);
    for ( int i = 0; i < 4; i++ ) {
        dataType->states[i] = String_clone(NUCLEOTIDES_STRING[i]);
    }
    dataType->encoding = _nucleotide_encoding;
    dataType->state    = _nucleotide_state;
    
    dataType->encoding_string = _nucleotide_encoding_string;
    dataType->state_string    = _state_string;
    
    dataType->state_count = _state_count;
    
    return dataType;
}

DataType *nucleotide_datatype(){
    return &SINGLETON_DATATYPE_NUCLEOTIDE;
}

int _nucleotide_encoding( DataType *datatype, char nuc){
    return NUCLEOTIDE_STATES[nuc];
}

int _nucleotide_encoding_string( DataType *datatype, const char *nuc){
    return NUCLEOTIDE_STATES[nuc[0]];
}

char _nucleotide_state( DataType *datatype, int encoding){
    return NUCLEOTIDES[encoding];
}

const char * _nucleotide_state_string( const DataType *datatype, int encoding){
    return NUCLEOTIDES_STRING[encoding];
}


#pragma mark-
#pragma mark Amino acid

DataType *new_AminoAcidDataType(){
    DataType *dataType = malloc(sizeof(DataType));
    assert(dataType);
    dataType->name = String_clone("Amino Acid");
    dataType->type = DATA_TYPE_AMINO_ACID;
    dataType->genetic_code = 0;
    dataType->stateCount = 20;
    dataType->symbolLength = 1;
    dataType->states = malloc(20*sizeof(char*));
    assert(dataType->states);
    for ( int i = 0; i < 20; i++ ) {
        dataType->states[i] = String_clone(AMINO_ACIDS_STRING[i]);
    }
    dataType->encoding = _aa_encoding;
    dataType->state    = _aa_state;
    
    dataType->encoding_string = _aa_encoding_string;
    dataType->state_string    = _state_string;
    
    dataType->state_count = _state_count;
    
    return dataType;
}

DataType *amino_acid_datatype(){
    return &SINGLETON_DATATYPE_AMINO_ACID;
}

int _aa_encoding( DataType *datatype, char aa){
    return AMINO_ACID_STATES[aa];
}

char _aa_state( DataType *datatype, int encoding){
    return AMINO_ACIDS[encoding];
}

int _aa_encoding_string( DataType *datatype, const char *aa){
    return AMINO_ACID_STATES[aa[0]];
}

const char * _aa_state_string( const DataType *datatype, int encoding){
    return AMINO_ACIDS_STRING[encoding];
}

#pragma mark-
#pragma mark Codon

DataType *new_CodonDataType( int genetic_code ){
    DataType *dataType = malloc(sizeof(DataType));
    assert(dataType);
    dataType->name = String_clone("Codon");
    dataType->type = DATA_TYPE_CODON;
    dataType->genetic_code = genetic_code;
    dataType->stateCount = NUMBER_OF_CODONS[genetic_code];
    dataType->symbolLength = 3;
    dataType->states = malloc(dataType->stateCount*sizeof(char*));
    assert(dataType->states);
    for ( int i = 0; i < dataType->stateCount; i++ ) {
        dataType->states[i] = String_clone(CODON_TRIPLETS[i]);
    }
    dataType->encoding = _codon_encoding;
    dataType->state    = _codon_state;
    
    dataType->encoding_string = _codon_encoding_string;
    dataType->state_string    = _codon_state_string;
    
    dataType->state_count = _state_count;
    
    return dataType;
}

DataType *codon_datatype( int genetic_code ){
    SINGLETON_DATATYPE_CODON.genetic_code = genetic_code;
    return &SINGLETON_DATATYPE_CODON;
}

int _codon_encoding( DataType *datatype, char nuc){
    return 0;
}

char _codon_state( DataType *datatype, int encoding){
    return ' ';
}

int _codon_state_count( DataType *datatype){
    return NUMBER_OF_CODONS[datatype->genetic_code];
}


void encoding_to_codon( int encoding, int genetic_code, char *codon ){
    if( encoding >= NUMBER_OF_CODONS[genetic_code] ){
        codon[0] = '-';
        codon[1] = '-';
        codon[2] = '-';
        return;
    }
    int count = 0;
    int i = 0;
    for ( ; i < 64; i++ ) {
        if( GENETIC_CODE_TABLES[genetic_code][i] == '*' ) continue;
        if(count == encoding ){
            break;
        }
        count++;
    }
    
    codon[0] = CODON_TRIPLETS[i][0];
    codon[1] = CODON_TRIPLETS[i][1];
    codon[2] = CODON_TRIPLETS[i][2];
}

int codon_to_encoding( const char *codon, int genetic_code ){
    DataType nuctype = SINGLETON_DATATYPE_NUCLEOTIDE;
    
    int n1 = nuctype.encoding(&nuctype, codon[0] );
    int n2 = nuctype.encoding(&nuctype, codon[1] );
    int n3 = nuctype.encoding(&nuctype, codon[2] );
    
    int encoding = 65;
    if ( n1 > 3 || n2 > 3 || n3 > 3 ) {
        return 65;
    }
    else {
        encoding = n1*16+n2*4+n3;
        
        for ( int j = 0; j < n1*16+n2*4+n3; j++ ) {
            if ( GENETIC_CODE_TABLES[genetic_code][j] == '*' ) {
                encoding--;
            }
        }
        
    }
    return encoding;
}

const char * _codon_state_string( const DataType *datatype, int encoding){
    if( encoding >= NUMBER_OF_CODONS[datatype->genetic_code] ){
        return "???";
    }
    int count = 0;
    int i = 0;
    for ( ; i < 64; i++ ) {
        if( GENETIC_CODE_TABLES[datatype->genetic_code][i] == '*' ) continue;
        if(count == encoding ){
            break;
        }
        count++;
    }
    return CODON_TRIPLETS[i];
}

int _codon_encoding_string( DataType *datatype, const char *codon){
    DataType nuctype = SINGLETON_DATATYPE_NUCLEOTIDE;
    
    int n1 = nuctype.encoding(&nuctype, codon[0] );
    int n2 = nuctype.encoding(&nuctype, codon[1] );
    int n3 = nuctype.encoding(&nuctype, codon[2] );
    
    const char *table = GENETIC_CODE_TABLES[datatype->genetic_code];
    int encoding = 65;
    if ( n1 > 3 || n2 > 3 || n3 > 3 ) {
        return 65;
    }
    else {
        encoding = n1*16+n2*4+n3;
        
        for ( int j = 0; j < n1*16+n2*4+n3; j++ ) {
            if ( table[j] == '*' ) {
                encoding--;
            }
        }
        
    }
    return encoding;
}
