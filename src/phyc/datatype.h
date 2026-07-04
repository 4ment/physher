// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_datatype_h
#define PhyC_datatype_h

#include "mjson.h"
#include "hashtable.h"

typedef enum datatype{ DATA_TYPE_NUCLEOTIDE, DATA_TYPE_CODON, DATA_TYPE_AMINO_ACID, DATA_TYPE_GENERIC} datatype;

static const double NUCLEOTIDE_AMBIGUITY_STATES[18][4] = {
	// A
    {1,0,0,0},
    // C
    {0,1,0,0},
    // G
    {0,0,1,0},
    // T
    {0,0,0,1},
    // U
    {0,0,0,1},
    // R
    {1,0,1,0},
    // Y
    {0,1,0,1},
    // M
    {1,1,0,0},
    // W
    {1,0,0,1},
    // S
    {0,1,1,0},
    // K
    {0,0,1,1},
    
    //"A", "C", "G", "T", "T", "AG", "CT", "AC", "AT", "CG", "GT",
	
    // B
    {0,1,1,1},
    // D
    {1,0,1,1},
    // H
    {1,1,0,1},
    // V
    {1,1,1,0},
    // N
    {1,1,1,1},
    // ?
    {1,1,1,1},
    // -
    {1,1,1,1},
    //"CGT", "AGT", "ACT", "ACG", "ACGT", "ACGT", "ACGT"
};


typedef struct DataType{
    char *name;
    datatype type;
    int stateCount;
    int symbolLength; //1 for nucletoide and aa. 3 for codon and whatever for general DataType
    char **states;
    Hashtable* ambiguities;
    int (*encoding)( const struct DataType *, char );
    char (*state)( const struct DataType *, int  );
    
    int (*encoding_string)( const struct DataType *, const char *);
    const char * (*state_string)( const struct DataType *, int  );
	void (*partial)(const struct DataType *, int, double*);
    
    int (*state_count)( struct DataType *  );
    int8_t genetic_code;
	int ref_count;
}DataType;

DataType* new_DataType_from_json(json_node* node, Hashtable* table);

DataType * clone_DataType(const DataType *dataType);

DataType *new_NucleotideDataType();

DataType *new_AminoAcidDataType();

DataType *new_CodonDataType(int genetic_code);

DataType *new_GenericDataType( const char* name, size_t count, const char **states) ;

void GenericDataType_add_ambiguity(DataType *dataType, const char* ambiguity, size_t count, const char** states);

void free_DataType( DataType *dataType);

//void encoding_to_codon( int encoding, int genetic_code, char *codon );

//int codon_to_encoding( const char *codon, int genetic_code );

DataType *nucleotide_datatype();

DataType *amino_acid_datatype();

DataType *codon_datatype( int genetic_code );


#endif
