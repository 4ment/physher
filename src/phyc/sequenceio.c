/*
 *  sequenceio.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/2/10.
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

#include "sequenceio.h"

#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "utils.h"
#include "mstring.h"
#include "sequence.h"
#include "filereader.h"

Sequences * readSequences( const char *infile ){
	int what = -1;
	StringBuffer *buffer  = new_StringBuffer(1000);
	FileReader *reader = new_FileReader(infile, 1000);
	
	while ( reader->read_line(reader) ) {
		if( String_is_empty(reader->line) ) continue;
		if ( String_start_with(reader->line, ">", false) ) {
			what = 0;
			break;
		}
		StringBuffer_assign_uppercase(reader->buffer, buffer);
		if ( String_i_start_with(buffer->c, "#NEXUS", true) ){
			what = 1;
			break;
		}
        
        char *ptr = reader->line;
        while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
        if( *ptr >= 48 && *ptr <=57  ){
            what = 2;
        }
	}
	free_FileReader(reader);
	free_StringBuffer(buffer);
	
	Sequences *seqs = NULL;
	switch ( what ) {
		case 0:
			seqs = readFasta(infile);
            break;
        case 1:
            seqs = readNexus(infile,0);
            break;
        case 2:
            seqs = readPhylip(infile);
            break;
			
		default:
            assert(0);
	}
	return seqs;
}

// Read FASTA file
Sequences * readFasta ( const char *infile ){
	Sequences  *seqs = new_Sequences(100);
		
	FileReader *reader = new_FileReader(infile, 1000);
	
	StringBuffer *seq  = new_StringBuffer(1000);
	StringBuffer *name = new_StringBuffer(100);
	
	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);
		if ( reader->buffer->length == 0) {
			continue;
		}
		
		if ( String_start_with(reader->line, ">", false) ) {
			if ( seq->length != 0) {
                StringBuffer_uppercase(seq);
				Sequences_add(seqs, new_Sequence( &name->c[1], seq->c));
				
				StringBuffer_empty(seq);
			}
			
			StringBuffer_set_string(name, reader->line);
		}
		else {
			StringBuffer_append_string(seq, reader->line );
		}
	}
	
	if ( seq->length != 0) {
        StringBuffer_uppercase(seq);
		Sequences_add(seqs, new_Sequence( &name->c[1], seq->c));
	}
	
	Sequences_pack(seqs);
	
	free_FileReader(reader);
	free_StringBuffer(name);
	free_StringBuffer(seq);
	
	return seqs;
}



// Read NEXUS file
// index: index of the dataset. There can be more than 1 BEGIN DATA block
Sequences * readNexus ( const char *infile, int index ){
	Sequences  *seqs = (Sequences *)malloc( sizeof(Sequences));
	assert(seqs);
	
	int nSeq = 0;
	int capacity = 10;
	seqs->seqs = (Sequence **)malloc( capacity * sizeof(Sequence *));
    assert(seqs->seqs);
	FileReader *reader = new_FileReader(infile, 1000);
	
	StringBuffer *buffer  = new_StringBuffer(1000);
	StringBuffer *bufferUpper  = new_StringBuffer(1000);
	StringBuffer *name = new_StringBuffer(100);
	
	const char *NEXUS_DATA_BLOCK = "BEGIN DATA";
	const char *NEXUS_CHARACTER_BLOCK = "BEGIN CHARACTER";
	
	int ntax = 0;
	int nchar = 0;
	bool matrixFlag = false;

	int counter = 0;
	
	while ( reader->read_line(reader) ) {
		StringBuffer_assign_uppercase( reader->buffer, bufferUpper );
		
		// BEGIN DATA block
		if ( String_i_start_with(bufferUpper->c, NEXUS_DATA_BLOCK, true) || String_i_start_with(bufferUpper->c, NEXUS_CHARACTER_BLOCK, true) ) {
            
			if( counter != index){				
				counter++;
				continue;
			}
			while ( reader->read_line(reader) ) {
				
				// MATRIX block
				if ( matrixFlag ) {
					if( String_is_empty(reader->line) ) continue;
					
					// Start of a comment
					if ( String_i_start_with(reader->line, "[", true) ){
						char *ptr = reader->line;
						while ( *ptr != ']') {
							if( *ptr == '\0' ) {
								break;
							}
							ptr++;
						}
						
						// multi-line comment
						if( *ptr != ']') {
							while ( reader->read_line(reader) ) {
								ptr = reader->line;
								while ( *ptr != ']') {
									if( *ptr == '\0' ) {
										break;
									}
									ptr++;
								}
								if( *ptr == ']') break;
							}
						}
						continue;
					}
					
					if ( String_i_start_with(reader->line, ";", true) ) {
						//matrixFlag = false;
						break;
					}
					
					
					char *ptr = reader->line;
                    while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
					if ( String_i_start_with(reader->line, "'", true) ){
						ptr++;
						while ( *ptr != '\'' ) {
							StringBuffer_append_char(name, *ptr);
							ptr++;
						}
                        ptr++;
					}
					else {
						while ( *ptr != ' ' && *ptr != '\t' ) {
							StringBuffer_append_char(name, *ptr);
							ptr++;
						}
					}
					
					while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
					while ( *ptr != '\0' ){
						if ( *ptr != ' ' && *ptr != '\t' ){
                            StringBuffer_append_char(buffer, toupper(*ptr));
                        }
						ptr++;
					}
					
					if( nSeq == capacity ){
						capacity *= 2;
						seqs->seqs = realloc( seqs->seqs, capacity * sizeof(Sequence *));
					}
					seqs->seqs[nSeq] = (Sequence *)malloc(sizeof(Sequence));
					assert(seqs->seqs[nSeq]);
					
					StringBuffer_trim(buffer);
					if ( String_end_with(buffer->c, ";", true)) {
						//matrixFlag = false;
						buffer->c[buffer->length-1] = '\0';
						buffer->length--;
						String_rtrim(buffer->c);
					}
                    
					StringBuffer_uppercase(buffer);
					seqs->seqs[nSeq]->seq = String_clone(buffer->c);
					seqs->seqs[nSeq]->name = String_clone(name->c);
					seqs->seqs[nSeq]->length = buffer->length;
                  		
                    if ( nchar != buffer->length ){
						fprintf(stderr, "nchar=%d but %s sequence length is %zu\n", nchar, name->c, buffer->length);
					}
					
					StringBuffer_empty(buffer);
					StringBuffer_empty(name);
					
					nSeq++;
					// for now we just exit when we have the sequences, we are only interested in them
					if ( matrixFlag == false ) break;
					
				}
				else {
					StringBuffer_assign_uppercase( reader->buffer, bufferUpper );
					
					if ( String_contains_str(bufferUpper->c, "NTAX") ){
						char *ptr = reader->line+String_index_of_str(bufferUpper->c, "NTAX")+4;
						while ( *ptr < 48 || *ptr > 57 ) ptr++;
						while ( *ptr >= 48 && *ptr <= 57 ){
							StringBuffer_append_char(buffer, *ptr);
							ptr++;
						}
						ntax = atoi(buffer->c);
						
//						capacity = ntax;
//						seqs->seqs = (Sequence **)malloc( capacity * sizeof(Sequence *));
//						assert(seqs->seqs);
						
						StringBuffer_empty(buffer);
					}
					if ( String_contains_str(bufferUpper->c, "NCHAR") ){
						char *ptr = reader->line+String_index_of_str(bufferUpper->c, "NCHAR")+5;
						while ( *ptr < 48 || *ptr > 57 ) ptr++;
						while ( *ptr >= 48 && *ptr <= 57 ){
							StringBuffer_append_char(buffer, *ptr);
							ptr++;
						}
						nchar = atoi(buffer->c);
						StringBuffer_empty(buffer);				
					}
					if( String_i_start_with(bufferUpper->c, "MATRIX", true) ){
						matrixFlag = true;
//						// ntax is missing
//						if ( seqs == NULL ) {
//							capacity = 100;
//							seqs->seqs = (Sequence **)malloc( capacity * sizeof(Sequence *));
//							assert(seqs->seqs);
//						}
					}
				}
				
			}
			break;
			
		}
		
	}
	
	
	free_FileReader(reader);
	free_StringBuffer(name);
	free_StringBuffer(buffer);
	free_StringBuffer(bufferUpper);
	
    if(nSeq == 0){
        error("Empty or unreadble nexus sequence file");
    }
    if( ntax!= 0 && ntax != nSeq){
        printf("ntax (%d) is not equal to number of sequences (%d) found in the nexus file\n", ntax, nSeq);
    }
    
	if ( nSeq != capacity ) seqs->seqs = realloc( seqs->seqs, nSeq * sizeof(Sequence *));
	seqs->size = nSeq;
	seqs->length = seqs->seqs[0]->length;
	
	return seqs;
}

// Read Phylip file
Sequences * readPhylip ( const char *infile ){
    Sequences  *seqs = new_Sequences(100);
    
    FileReader *reader = new_FileReader(infile, 1000);
    
    StringBuffer *seq  = new_StringBuffer(1000);
    StringBuffer *name = new_StringBuffer(100);
    
    int ntax = 0;
    int nchar = 0;
    bool first = true;
    int i ;
    
    int count = 0;
    int index = 0;
    
    while ( reader->read_line(reader) ) {
        StringBuffer_trim(reader->buffer);
        if ( reader->buffer->length == 0) {
            continue;
        }
        char *ptr = reader->line;
        
        // read first line and get ntax and nchar
        if( first ){
            while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
            while ( *ptr >= 48 && *ptr <=57 ){
                StringBuffer_append_char(name, *ptr);
                ptr++;
            }
            ntax = atoi(name->c);
            
            StringBuffer_empty(name);
            while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
            while ( *ptr >= 48 && *ptr <=57 ){
                StringBuffer_append_char(name, *ptr);
                ptr++;
            }
            nchar = atoi(name->c);
            first = false;
        }
        else if( count == ntax ){
            StringBuffer_empty(seq);
            while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
            while ( *ptr != ' ' && *ptr != '\t' ){
                StringBuffer_append_char(seq, *ptr);
                ptr++;
            }
            StringBuffer_uppercase(seq);
            
            seqs->seqs[index]->seq = String_append_string(seqs->seqs[index]->seq, seq->c);
            seqs->seqs[index]->length = strlen(seqs->seqs[index]->seq);
            index++;
            if( index == ntax ) index = 0;
        }
        else {
            while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
            
            StringBuffer_empty(name);
            while ( *ptr != ' ' && *ptr != '\t' ){
                StringBuffer_append_char(name, *ptr);
                ptr++;
            }
            
            StringBuffer_empty(seq);
            while ( *ptr == ' ' || *ptr == '\t' ) ptr++;
            while ( *ptr != ' ' && *ptr != '\t' ){
                StringBuffer_append_char(seq, *ptr);
                ptr++;
            }
            
            StringBuffer_uppercase(seq);
            Sequences_add(seqs, new_Sequence( name->c, seq->c));
            count++;
        }
    }
    
    Sequences_pack(seqs);
    
    free_FileReader(reader);
    free_StringBuffer(name);
    free_StringBuffer(seq);
    
    for ( i = 0; i < seqs->size; i++ ) {
        if( nchar != seqs->seqs[i]->length ){
            fprintf(stderr,"nchar (%d) is not equal to sequence length (%d) in the phylip file\n", nchar, seqs->seqs[i]->length);
            exit(1);
        }
    }
    seqs->size = count;
    seqs->length = seqs->seqs[0]->length;
    
    return seqs;
}

void Sequences_save_fasta( const Sequences *sequences, const char *filename ){
    FILE *pfile = NULL;
    if( filename == NULL ){
        pfile = stdout;
    }
    else {
        pfile = fopen(filename, "w");
    }
	
	for ( int i = 0; i < sequences->size; ++i ) {
		fprintf(pfile, ">%s\n", sequences->seqs[i]->name);
		fprintf(pfile, "%s\n", sequences->seqs[i]->seq);
	}
	fclose(pfile);
}

void Sequences_save_nexus( const Sequences *sequences, const char *filename ){
	Sequences_save_nexus_with_comment(sequences, filename, NULL);
}

void Sequences_save_nexus_with_comment( const Sequences *sequences, const char *filename, char *comment ){
	FILE *pfile = fopen(filename, "w");
	fprintf(pfile, "#NEXUS\n\n\n");
    if( comment != NULL )  fprintf(pfile, "%s\n\n\n",comment);
	fprintf(pfile, "Begin DATA;\n");
	fprintf(pfile, "\tDimensions ntax=%d nchar=%d;\n", sequences->size, sequences->length);
	fprintf(pfile, "\tFormat datatype=%s gap=-;\n", sequences->datatype->name);
	fprintf(pfile, "\tMatrix\n");
    
	int max_name = 0;
    bool containsSpace = false;
	for ( int i = 0; i < sequences->size; ++i ) {
		max_name = imax(max_name, strlen(sequences->seqs[i]->name));
        for (int j = 0; j < strlen(sequences->seqs[i]->name); j++) {
            if(sequences->seqs[i]->name[j] == ' '){
                containsSpace = true;
                break;
            }
        }
	}
    
	for ( int i = 0; i < sequences->size; ++i ) {
        if(containsSpace){
            fprintf(pfile, "\'%s\'", sequences->seqs[i]->name);
        }
        else{
            fprintf(pfile, "%s", sequences->seqs[i]->name);
        }
		for ( int s = 0; s < (max_name-strlen(sequences->seqs[i]->name)+6); s++ ) {
			fprintf(pfile, "%c", ' ');
		}
		fprintf(pfile, "%s\n", sequences->seqs[i]->seq);
	}
	fprintf(pfile, "\n;\nEnd;\n");
	fclose(pfile);
}

void Sequences_save_phylip( const Sequences *sequences, const char *filename ){
	FILE *pfile = fopen(filename, "w");
    fprintf(pfile, " %d %d\n", sequences->size,sequences->length);
	for ( int i = 0; i < sequences->size; ++i ) {
		fprintf(pfile, "%s  ", sequences->seqs[i]->name);
		fprintf(pfile, "%s\n", sequences->seqs[i]->seq);
	}
	fclose(pfile);
}

Sequences* new_Sequences_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"file",
		"sequences"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	json_node* file_node = get_json_node(node, "file");
	json_node* sequences_node = get_json_node(node, "sequences");
	Sequences* sequences = NULL;
	
	if(sequences_node != NULL){
		sequences = new_Sequences(sequences_node->child_count);
		for (int i = 0; i < sequences_node->child_count; i++) {
			json_node* child = sequences_node->children[i];
			char* taxon = child->key;
			char* sequence = child->value;
			Sequences_add(sequences, new_Sequence(taxon, sequence));
		}
	}
	else if(file_node != NULL){
		const char* filename = (const char*)file_node->value;
		sequences = readSequences(filename);
	}
	else{
		fprintf(stderr, "No `sequences' input (use file option)\n");
		exit(1);
	}
	return sequences;
}
