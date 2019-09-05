/*
 *  nexus.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/1/10.
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

#include "treeio.h"

#include <string.h>
#include <strings.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>

#include "utils.h"
#include "mstring.h"
#include "tree.h"
#include "filereader.h"
#include "modelselection.h"

#include <unistd.h> // for sleep


const char *TREE_TAG = "TREE";
const char *NEXUS_TREE_BLOCK = "BEGIN TREES";
const char *NEXUS_TREE_TRANSLATE_BLOCK = "TRANSLATE";

void print_tree_extended( FILE *pf, const Node *n, char **info ){
	if( n == NULL ) return;
	if( n->left != NULL ) fprintf(pf, "(");
	else{
		fprintf(pf, "%s[&%s]:%f", n->name, info[n->postorder_idx], n->distance->value);
		return;
	}
	print_tree_extended( pf, n->left, info );
	fprintf(pf, ",");
	print_tree_extended( pf, n->right, info );
	fprintf(pf, ")");
	if( n->parent != NULL ) fprintf(pf, "[&%s]:%f", info[n->postorder_idx], n->distance->value);
	//else fprintf(pf, "\n");	
}


void Tree_print_nexus_header_figtree( FILE *pf, Tree *tree ){
	Tree_print_nexus_header_figtree_Taxa(pf, tree);
	Tree_print_nexus_header_figtree_BeginTrees(pf, tree);
}

// Taxlabels are in preorder
void Tree_print_nexus_taxa_block( FILE *pf, Tree *tree ){
    
    fprintf(pf, "Begin taxa;\n");
    fprintf(pf, "\tDimensions ntax=%d;\n", Tree_tip_count(tree) );
    fprintf(pf, "\tTaxlabels\n");
    
    Node **nodes  = get_tips( tree, PREORDER );
    
    for ( int i = 0; i < Tree_tip_count(tree); i++) {
        if( Node_isleaf(nodes[i]) ) fprintf(pf, "\t\t%s\n", nodes[i]->name);
    }
    fprintf(pf, "\t\t;\nEnd;\n\n");
    free(nodes);
}

void Tree_print_nexus_header_figtree_Taxa( FILE *pf, Tree *tree ){
    
    fprintf(pf, "#NEXUS\n\nBegin taxa;\n");
    fprintf(pf, "\tDimensions ntax=%d;\n", Tree_tip_count(tree) );
    fprintf(pf, "\tTaxlabels\n");
    
    Node **nodes  = get_tips( tree, PREORDER );
    
    for ( int i = 0; i < Tree_tip_count(tree); i++) {
        if( Node_isleaf(nodes[i]) ) fprintf(pf, "\t\t%s\n", nodes[i]->name);
    }
    fprintf(pf, "\t\t;\nEnd;\n\n");
    free(nodes);
}

// Translation indexes are in preorder
void Tree_print_nexus_header_figtree_BeginTrees( FILE *pf, Tree *tree ){
	Node **nodes  = Tree_get_nodes(tree, PREORDER);
	fprintf(pf, "Begin trees;\n");
	fprintf(pf, "\tTranslate\n");
    
    int count = 1;
    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
        if( Node_isleaf(nodes[i]) ){
            fprintf(pf, "\t\t%d %s", count, Node_name(nodes[i]));
			if( count != Tree_tip_count(tree) ) fprintf(pf, ",");
			fprintf(pf, "\n");
            count++;
        }
    }
	
	fprintf(pf, "\t\t;\n");
}

static void _Tree_print_nexus_with_annotation_aux( FILE *pf, Tree *tree, const Node *n, int *count, bool time ){
    if( n == NULL ) return;
	if( !Node_isleaf(n) ) fprintf(pf, "(");
	else {
        if( n->annotation != NULL && Hashtable_length(n->annotation) > 0 ){
            StringBuffer *buff = new_StringBuffer(10);
            Hashtable_init_iterator(n->annotation);
            HashEntry *entry = NULL;
            while ( (entry = Hashtable_next(n->annotation) ) != NULL ) {
                char *value = (char*)HashEntry_value(entry);
                char *key   = (char*)HashEntry_key(entry);
                StringBuffer_append_strings(buff, 4, key,"=",value, ",");
            }
            StringBuffer_chop(buff);
            if(time){
                fprintf(pf, "%d:[&%s]%.20f", ++(*count), buff->c, (n->parent->height->value - n->height->value) );
            }
            else {
                fprintf(pf, "%d:[&%s]%.20f", ++(*count), buff->c, n->distance->value );
            }
            //fflush(pf);
            free_StringBuffer(buff);
        }
        else {
            if(time){
                fprintf(pf, "%d:%.20f", ++(*count), (n->parent->height->value - n->height->value) );
            }
            else{
                fprintf(pf, "%d:%.20f", ++(*count), n->distance->value );
            }
        }
		return;
	}
	
	_Tree_print_nexus_with_annotation_aux( pf, tree, n->left, count, time );
	
	fprintf(pf, ",");
	
	_Tree_print_nexus_with_annotation_aux( pf, tree, n->right, count, time );
    
    if( n->annotation != NULL && Hashtable_length(n->annotation) > 0 ){
        StringBuffer *buff = new_StringBuffer(10);
        Hashtable_init_iterator(n->annotation);
        HashEntry *entry = NULL;
        while ( (entry = Hashtable_next(n->annotation) ) != NULL ) {
            char *value = (char*)HashEntry_value(entry);
            char *key   = (char*)HashEntry_key(entry);
            StringBuffer_append_strings(buff, 4, key,"=",value, ",");
        }
        StringBuffer_chop(buff);
        if( Node_isroot(n) ){
            fprintf(pf, ")[&%s];", buff->c );
        }
        else {
            if(time){
                fprintf(pf, "):[&%s]%f", buff->c, Node_time_elapsed((Node*)n) );
            }
            else{
                fprintf(pf, "):[&%s]%f", buff->c, Node_distance((Node*)n) );
            }
        }
        free_StringBuffer(buff);
    }
    else {
        if( Node_isroot(n) ) fprintf(pf, ");" );
        else{
            if(time){
                fprintf(pf, "):%f", Node_time_elapsed((Node*)n) );
            }
            else {
                fprintf(pf, "):%f", Node_distance((Node*)n) );
            }
        }
    }
}

void Tree_print_nexus_with_annotation( FILE *pf, Tree *tree ){
    int cunt = 0;
    _Tree_print_nexus_with_annotation_aux(pf, tree, Tree_root(tree), &cunt, true);
}

void Tree_print_nexus_with_annotation2( FILE *pf, Tree *tree, bool time ){
    int cunt = 0;
    _Tree_print_nexus_with_annotation_aux(pf, tree, Tree_root(tree), &cunt, time);
}

static void _Tree_print_nexus_aux( FILE *pf, Tree *tree, const Node *n, int *count ){
    if( n == NULL ) return;
	if( !Node_isleaf(n) ) fprintf(pf, "(");
	else {
		if(Tree_is_time_mode(tree)){
        	fprintf(pf, "%d:%f", ++(*count), (n->parent->height->value - n->height->value) );
		}
		else{
			fprintf(pf, "%d:%f", ++(*count), Node_distance(n) );
		}
		return;
	}
	
	_Tree_print_nexus_aux( pf, tree, n->left, count );
	
	fprintf(pf, ",");
	
	_Tree_print_nexus_aux( pf, tree, n->right, count );
	
    if( Node_isroot(n) ){
        fprintf(pf, ");" );
    }
    else {
		if(Tree_is_time_mode(tree)){
        	fprintf(pf, "):%f", Node_distance(n) );
		}
		else{
			fprintf(pf, "):%f", Node_time_elapsed((Node*)n) );
		}
    }
}

void Tree_print_nexus( FILE *pf, Tree *tree ){
    int cunt = 0;
    _Tree_print_nexus_aux(pf, tree, Tree_root(tree), &cunt);
}



void Tree_print_newick_subtree( FILE *pf, bool time, const Node *n, bool internal ){
	if( n == NULL ) return;
	if( n->left != NULL ) fprintf(pf, "(");
	else{
		if(time){
			fprintf(pf, "%s:%.8f", n->name, Node_time_elapsed((Node*)n) );
		}
		else{
			fprintf(pf, "%s:%.8f", n->name, n->distance->value);
		}
		return;
	}
	Tree_print_newick_subtree( pf, time, n->left, internal );
	fprintf(pf, ",");
	Tree_print_newick_subtree( pf, time, n->right, internal );
	fprintf(pf, ")");
	if( n->parent != NULL ){
		if( internal ) fprintf(pf, "%s", n->name);
		if(time){
			fprintf(pf, ":%.8f", Node_time_elapsed((Node*)n));
		}
		else{
			fprintf(pf, ":%.8f", n->distance->value);
		}
	}
	else fprintf(pf, ";");
}

void Tree_print_newick( FILE *pf, Tree *tree, bool internal ){
    Tree_print_newick_subtree(pf, Tree_is_time_mode(tree), Tree_root(tree), internal);
}

static void _Tree_print_height_newick_aux( FILE *pf, Tree *tree, const Node *n, bool internal ){
	if( n == NULL ) return;
	if( !Node_isleaf(n) ) fprintf(pf, "(");
	else {
		fprintf(pf, "%s:%.8f", n->name, Node_time_elapsed((Node*)n) );
		return;
	}
	
	_Tree_print_height_newick_aux( pf, tree, n->left, internal );
	fprintf(pf, ",");
	_Tree_print_height_newick_aux( pf, tree, n->right, internal );
	if( Node_isroot(n) ){
		fprintf(pf, ");" );
	}
	else {
		fprintf(pf, ")%s:%f", (internal?n->name:""), Node_time_elapsed((Node*)n) );
	}
}

void Tree_print_height_newick( FILE *pf, Tree *tree, bool internal ){
	_Tree_print_height_newick_aux(pf, tree, Tree_root(tree), internal);
}


#pragma mark -
// MARK: read tree

char *readTree( const char *infile ){
    
    enum FORMAT { FASTA, NEXUS, TREE, PHYLIP };
    enum FORMAT format = NEXUS;
    
    FileReader *reader = new_FileReader(infile, 1000);
    
    while ( reader->read_line(reader) ) {
        if( String_start_with(reader->line, "#NEXUS", false) ) format = NEXUS;
        else if( String_start_with(reader->line, ">", false) ) format = FASTA;
        else if( String_start_with(reader->line, "(", true) ) format = TREE;
        else {
            error("Incompatible file format (Only FASTA, NEXUS format)\n");
        }
        break;
    }
    
    free_FileReader(reader);
    switch(format){
        case NEXUS:{
            return readNexusFirstTree(infile);
        }
        case FASTA:
            return readNewickTree(infile);
            break;
        case TREE:
            return readNewickTree(infile);
            break;
        default: return NULL;
    }
}

char **readTrees( const char *infile, int *trees_count ){
    
    return readTreesRange(infile, trees_count, 0, INT_MAX);
}

char **readTreesRange( const char *infile, int *trees_count, int start, int end ){
    
    enum FORMAT { FASTA, NEXUS, TREE, PHYLIP };
    enum FORMAT format = NEXUS;
    
    FileReader *reader = new_FileReader(infile, 1000);
    
    while ( reader->read_line(reader) ) {
        if( String_start_with(reader->line, "#NEXUS", false) ) format = NEXUS;
        else if( String_start_with(reader->line, "(", true) ) format = TREE;
        else {
            error("Incompatible file format (Only FASTA, NEXUS format)\n");
        }
        break;
    }
    
    free_FileReader(reader);
    switch(format){
        case NEXUS:{
            return readNexusTreesRange(infile, trees_count, start, end);
        }
        case TREE:
            return readNewickTrees(infile, trees_count);
            break;
        default: return NULL;
    }
}

char * readNexusFirstTree( const char *infile ){
	int count = 0;
	char **trees = readNexusTrees( infile, &count, NULL, 1 );
	if ( count == 1 ) {
		char *tree = trees[0];
		free(trees); trees = NULL;
		return tree;
	}
	else {
		return NULL;
	}
}

char * readNexusTreeAt( const char *infile, const int index ){
	int count = 0;
	char str[1000];
	sprintf(str,"%d-%d", index, index);
	char **trees = readNexusTrees( infile, &count, str, 1 );
	if ( count == 1 ) {
		char *tree = trees[0];
		free(trees); trees = NULL;
		return tree;
	}
	else {
		return NULL;
	}
}

char ** readNexusRangeTree( const char *infile, int *count, char *range ){
	return readNexusTrees( infile, count, ",", 100000000 );
}

char ** readNexusNFirstTree( const char *infile, int *count, int max ){
	return readNexusTrees( infile, count, NULL, max );
}

char ** parseNexusTranslateBlock( FileReader *reader, StringBuffer *buffer, unsigned *count ){
    unsigned capcacity = 50;
    char **names = (char**)calloc( capcacity, sizeof(char*) );
    assert(names);
    
    bool end_block = false;
    
    while ( reader->read_line(reader) ) {
        if ( String_equals(reader->line, ";", true)) {
            break;
        }
        // e.g. 2 A/Singapore/2012[,;]
        else {
            char * line = reader->line;
            
            line = String_trim(line);
            StringBuffer_empty(buffer);
            
            
            while ( !isspace(*line ) ) {
                StringBuffer_append_char(buffer, *line);
                line++;
            }
            
            
            while ( isspace(*line ) ) {line++;}
            
            if ( *count == capcacity-1) {
                capcacity *= 2;
                names = realloc(names, capcacity*sizeof(char*) );
            }
            
            if ( String_end_with(line, ";", false)) {// already trimmed
                line[strlen(line)-1] = '\0';
                while ( isspace(line[strlen(line)-1] ) ) {
                    line[strlen(line)-1] = '\0';
                }
                end_block = true;
            }
            else if ( String_end_with(line, ",", false)) {
                line[strlen(line)-1] = '\0';
                while ( isspace(line[strlen(line)-1] ) ) {
                    line[strlen(line)-1] = '\0';
                }
            }
            
            //fprintf(stderr, "%d %s ",*count, line);
            // quotes
            if ( line[ strlen(line)-1 ] == '\'' || line[ strlen(line)-1 ] == '"' ) line[ strlen(line)-1 ] = '\0';
            if ( line[0] == '\'' || line[0] == '"' ) line++;
            //fprintf(stderr, "%s\n",line);
            names[*count] = String_clone(line);
            (*count)++;
            
            if ( end_block ) {
                break;
            }
            
        }
        
    }
    assert(*count > 0);
    names = realloc(names, *count*sizeof(char*) );
    return names;
}

Hashtable * parseNexusTranslateBlock2( FileReader *reader){
    Hashtable *hash = new_Hashtable_string(10);
    StringBuffer *buffer = new_StringBuffer(10);
    
    bool end_block = false;
    
    while ( reader->read_line(reader) ) {
        String_trim(reader->line);
        if ( reader->line[0] == ';' ) {
            break;
        }
        // e.g. 2 A/Singapore/2012[,;]
        else {
            char * line = reader->line;
            
            while ( isspace(*line ) ) {line++;}
            
            StringBuffer_empty(buffer);
            while ( !isspace(*line ) ) {
                StringBuffer_append_char(buffer, *line);
                line++;
            }
            
            char *key = StringBuffer_tochar(buffer);
            
            while ( isspace(*line ) ) {line++;}
            
            StringBuffer_empty(buffer);
            while ( *line != ';' && *line != ',' && *line != '\0' ) {
                StringBuffer_append_char(buffer, *line);
                line++;
            }
            if( *line == ';'){
                end_block = true;
            }
            
            StringBuffer_trim(buffer);
            if ( buffer->c[0] == '\'' || buffer->c[0] == '"' ){
                StringBuffer_chop(buffer);
                StringBuffer_cut(buffer);
            }
            //printf("%s | %s\n", key, buffer->c);
            Hashtable_add(hash, key, String_clone(buffer->c));
            
            if ( end_block ) {
                break;
            }
            
        }
        
    }
    free_StringBuffer(buffer);
    return hash;
}

char * insertNexusNames( char *tree, char ** const names ){
    StringBuffer *buffer = new_StringBuffer(strlen(tree)+100);
    StringBuffer *ibuffer = new_StringBuffer(20);
    char *ptr = tree;
    char *previous = tree;
    
    StringBuffer_append_char(buffer, *ptr); // it is never a taxon
    ++ptr;
    while ( *ptr != '\0' ) {
        if( *ptr == '[' ) {
            StringBuffer_append_char(buffer, *ptr);
            while ( *ptr != ']') {
                ptr++;
                previous++;
                StringBuffer_append_char(buffer, *ptr);
            }
            previous++;
            ptr++;
            continue;
        }
        else
            if ( (*previous == '(' || *previous == ',') && (*ptr >= 48  && *ptr <= 57) ) {
                StringBuffer_empty(ibuffer);
                while ( *ptr >= 48  && *ptr <= 57 ) {
                    StringBuffer_append_char(ibuffer, *ptr);
                    ptr++; previous++;
                }
                int index = atoi(ibuffer->c);
                StringBuffer_append_string(buffer, names[index-1]);
            }
        
        StringBuffer_append_char(buffer, *ptr);
        
        ptr++; previous++;
    }
    free_StringBuffer(ibuffer);
    ptr = NULL;
    previous = NULL;
    tree = realloc(tree, (buffer->length+1)*sizeof(char));
    strcpy(tree, buffer->c);
    //tree = StringBuffer_assign(buffer, tree);
    free_StringBuffer(buffer);
    return tree;
}

char * insertNexusNames2( char *tree, Hashtable *hash ){
    StringBuffer *buffer = new_StringBuffer(strlen(tree)+100);
    StringBuffer *ibuffer = new_StringBuffer(20);
    char *ptr = tree;
    char *previous = tree;
    
    StringBuffer_append_char(buffer, *ptr); // it is never a taxon
    ++ptr;
    while ( *ptr != '\0' ) {
        if( *ptr == '[' ) {
            StringBuffer_append_char(buffer, *ptr);
            while ( *ptr != ']') {
                ptr++;
                previous++;
                StringBuffer_append_char(buffer, *ptr);
            }
            previous++;
            ptr++;
            continue;
        }
        else {
            if ( (*previous == '(' || *previous == ',') && *ptr != '(' ) {
                StringBuffer_empty(ibuffer);
                while ( *ptr != ':' ) {
                    StringBuffer_append_char(ibuffer, *ptr);
                    ptr++; previous++;
                }
                char *taxon = Hashtable_get(hash, ibuffer->c);
                //printf("taxon %s - %s\n", ibuffer->c, taxon);
                StringBuffer_append_string(buffer, taxon);
            }
        }
        StringBuffer_append_char(buffer, *ptr);
        
        ptr++; previous++;
    }
    
    char *t = String_clone(buffer->c);
    free_StringBuffer(ibuffer);
    free_StringBuffer(buffer);

    return t;
}


void insertNexusNamesStringBuffer( StringBuffer *tree, char ** const names ){
	StringBuffer *buffer = new_StringBuffer(tree->length+100);
	StringBuffer *ibuffer = new_StringBuffer(20);
	char *ptr = tree->c;
	char *previous = tree->c;
	
	StringBuffer_append_char(buffer, *ptr); // it is never a taxon
	++ptr;
	while ( *ptr != '\0' ) {
		if ( (*previous == '(' || *previous == ',') && (*ptr >= 48  && *ptr <= 57) ) {
			StringBuffer_empty(ibuffer);
			while ( *ptr >= 48  && *ptr <= 57 ) {
				StringBuffer_append_char(ibuffer, *ptr);
				ptr++; previous++;
			}
			int index = atoi(ibuffer->c);
			StringBuffer_append_string(buffer, names[index-1]);
		}
		
		StringBuffer_append_char(buffer, *ptr);
		
		ptr++; previous++;
	}
	free_StringBuffer(ibuffer);
	
	// exchange pointers
	ptr = buffer->c;
	buffer->c = tree->c;
	tree->c = ptr;
	
	free_StringBuffer(buffer);
}

TreeFileIterator *new_TreeFileIterator( const char *filename ){
    FileReader *reader = new_FileReader(filename, 100);
    TreeFileIterator *iter = malloc(sizeof(TreeFileIterator));
    assert(iter);
    
    iter->reader = reader;
    iter->translation = NULL;
	iter->more = false;
    StringBuffer *buffer = new_StringBuffer(200);
    
    if( (iter->more = reader->read_line(iter->reader)) ){
        if(strcasecmp(reader->line, "#nexus") == 0){
            while ( reader->read_line(iter->reader) ) {
                // PARSE TREE BLOCK
                if ( String_i_start_with(reader->line, NEXUS_TREE_BLOCK, true) ) {
                    // PARSE TRANSLATE BLOCK
                    while ( reader->read_line(reader) ) {
                        if ( String_i_start_with(reader->line, NEXUS_TREE_TRANSLATE_BLOCK, true) ) {
                            iter->translation = parseNexusTranslateBlock2(reader);
                        }
                        else if ( String_i_start_with(reader->line, TREE_TAG, true) ) break;
                    }
                    break;
                }
            }
        }
		iter->more = true;
    }
    free_StringBuffer(buffer);
    return iter;
}

void free_TreeFileIterator(TreeFileIterator *iter ){
    
    free_FileReader(iter->reader);
    if( iter->translation != NULL){
        free_Hashtable(iter->translation);
    }
    free(iter);
}

char * TreeFileIterator_next_tree( TreeFileIterator *iter ){
    char *t = NULL;
    if( iter->more ){
        String_trim(iter->reader->line);
		if(strlen(iter->reader->line) > 0 && iter->reader->line[0] == '('){
			t = String_clone(iter->reader->line);
			iter->more = iter->reader->read_line(iter->reader);
			if(iter->more && strlen(iter->reader->line) == 0) iter->more = false;
		}
        else if( !String_i_start_with(iter->reader->line, "end;", true) ){
            char *ptr = iter->reader->line;

            // check comments
            // not handling multiline comments
            bool comment = false;
            while ( *ptr != '(' || comment ) {
                if ( *ptr == '[' ) {
                    comment = true;
                }
                else if( comment && *ptr == ']' ) {
                    comment = false;
                }
                ptr++;
            }
            
            // remove comments and check if there is ; at the end
            char *pch = &ptr[ strlen(ptr)-1 ];
            while( *pch != ';' && *pch != *ptr ){
                pch--;
            }
            
            if ( *pch == *ptr ) {
                error("missing ; in tree file\n");
            }
            *pch = '\0';
            
            //print_Hashtable(stdout, iter->translation);
            t = insertNexusNames2(ptr, iter->translation);
            
            // get the next tree if any
            while ( iter->reader->read_line(iter->reader) ) {
				if ( String_i_start_with(iter->reader->line, TREE_TAG, true) ){
					iter->more = true;
					break;
				}
            }
        }
    }
    
    return t;
}

char ** readNexusTreesRange( const char *infile, int *trees_count, int start, int end ){
    
    int nRange = 0;
    char **positions = NULL;
    
    const char *TREE_TAG = "TREE";
    const char *NEXUS_TREE_BLOCK = "BEGIN TREES";
    const char *NEXUS_TREE_TRANSLATE_BLOCK = "TRANSLATE";
    
    //const char *NEXUS_DATA_BLOCK = "BEGIN DATA";
    
    int trees_capcacity = 100;
    char **trees = (char**)calloc( trees_capcacity, sizeof(char*) );
    assert(trees);
    
    unsigned names_count = 0;
    char **names = NULL;
    
    FileReader *reader = new_FileReader(infile, 1000);
    StringBuffer *buffer = new_StringBuffer(200);
    
    *trees_count = 0;
    int trees_index = 0;
    //printf("start %d %d\n", start, end);
    while ( reader->read_line(reader) ) {
        
        // PARSE TREE BLOCK
        if ( String_i_start_with(reader->line, NEXUS_TREE_BLOCK, true) ) {
            while ( reader->read_line(reader) ) {
                // PARSE TRANSLATE BLOCK
                if ( String_i_start_with(reader->line, NEXUS_TREE_TRANSLATE_BLOCK, true) ) {
                    names = parseNexusTranslateBlock(reader, buffer, &names_count);
                }
                
                if ( String_i_start_with(reader->line, TREE_TAG, true) ) {
                    //printf("== %s\n",reader->line);
                    if( trees_index >= start){
                        if( trees_index > end) break;
                        
                        // check comments
                        // not handling multiline comments
                        char *ptr = reader->line;
                        bool comment = false;
                        while ( *ptr != '(' || comment ) {
                            if ( *ptr == '[' ) {
                                comment = true;
                            }
                            else if( comment && *ptr == ']' ) {
                                comment = false;
                            }
                            ptr++;
                        }
                        
                        // remove comments and check if there is ; at the end
                        char *pch = &ptr[ strlen(ptr)-1 ];
                        while( *pch != ';' && *pch != *ptr ){
                            pch--;
                        }
                        
                        if ( *pch == *ptr ) {
                            error("missing ; in tree file\n");
                        }
                        *pch = '\0';
                        
                        if ( *trees_count == trees_capcacity) {
                            trees_capcacity *= 2;
                            trees = realloc(trees, trees_capcacity*sizeof(char*) );
                        }
                        
                        trees[*trees_count] = String_clone(ptr);
                        ++(*trees_count);
                    }
                    trees_index++;
                }
                else if( String_i_start_with(reader->line, "END", true)){
                    break;
                }
            }
        }
    }
    
    free_FileReader(reader);
    free_StringBuffer(buffer);
    
    //fprintf(stderr, "count %d\n", *trees_count);
    
    if ( names_count != 0 ) {
        for (int i = 0; i < *trees_count; i++) {
            trees[i] = insertNexusNames(trees[i], names);
        }
        for (int i = 0; i < names_count; i++) {
            free(names[i]); names[i] = NULL;
        }
        free(names); names = NULL;
    }
    assert(*trees_count>0);
    if( trees_capcacity != *trees_count){
        trees = realloc(trees, *trees_count*sizeof(char*) );
    }
    
    if ( positions != NULL){
        for ( int i = 0; i < nRange; i++) {
            free(positions[i]);
        }
        free(positions);
    }
    
    if ( *trees_count == 0 ) {
        free(trees);
        trees = NULL;
    }
    return trees;
}

char ** readNexusTrees( const char *infile, int *trees_count, char *range, int max ){
    
    int nRange = 0;
    char **positions = NULL;
    
    int start = -1;
    int end   = -1;
    
    if ( range != NULL) {
        positions = String_split_char(range, ',', &nRange);
        if ( positions == NULL) {
            positions = String_split_char(range, '-', &nRange);
            start = atoi( positions[0] );
            end   = atoi( positions[1] );
        }
    }
    
    int tree_index = 0; // index in the nexus file of the current tree
    
    const char *TREE_TAG = "TREE";
    const char *NEXUS_TREE_BLOCK = "BEGIN TREES";
    const char *NEXUS_TREE_TRANSLATE_BLOCK = "TRANSLATE";
    
    //const char *NEXUS_DATA_BLOCK = "BEGIN DATA";
    
    int trees_capcacity = imin( (end-start+1), max);
    char **trees = (char**)calloc( trees_capcacity, sizeof(char*) );
    assert(trees);
    
    unsigned names_count = 0;
    char **names = NULL;
    
    FileReader *reader = new_FileReader(infile, 1000);
    StringBuffer *buffer = new_StringBuffer(200);
    
    while ( reader->read_line(reader) ) {
        
        // PARSE TREE BLOCK
        if ( String_i_start_with(reader->line, NEXUS_TREE_BLOCK, true) ) {
            while ( reader->read_line(reader) ) {
                // PARSE TRANSLATE BLOCK
                if ( String_i_start_with(reader->line, NEXUS_TREE_TRANSLATE_BLOCK, true) ) {
                    names = parseNexusTranslateBlock(reader, buffer, &names_count);
                }
                
                if ( String_i_start_with(reader->line, TREE_TAG, true) ) {
                    if ( range != NULL){
                        if ( start == -1 ) {
                            int i = 0;
                            for ( ; i < nRange; i++ ) {
                                int p = atoi(positions[i]);
                                if ( p == tree_index ) {
                                    break;
                                }
                            }
                            
                            if ( i == nRange ) {
                                tree_index++;
                                continue;
                            }
                        }
                        else {
                            if ( !( tree_index >= start   && tree_index <= end) ){
                                tree_index++;
                                continue;
                            }
                        }
                    }
                    
                    // check comments
                    // not handling multiline comments
                    char *ptr = reader->line;
                    bool comment = false;
                    while ( *ptr != '(' || comment ) {
                        if ( *ptr == '[' ) {
                            comment = true;
                        }
                        else if( comment && *ptr == ']' ) {
                            comment = false;
                        }
                        ptr++;
                    }
                    
                    // remove comments and check if there is ; at the end
                    char *pch = &ptr[ strlen(ptr)-1 ];
                    while( *pch != ';' && *pch != *ptr ){
                        pch--;
                    }
                    
                    if ( *pch == *ptr ) {
                        error("missing ; in tree file\n");
                    }
                    *pch = '\0';
                    
                    if ( *trees_count == trees_capcacity) {
                        trees_capcacity *= 2;
                        trees = realloc(trees, trees_capcacity*sizeof(char*) );
                    }
                    
                    trees[*trees_count] = String_clone(ptr);
                    
                    //fprintf(stderr, "\n%s\n", ptr);
                    
                    ++(*trees_count);
                    
                    if ( *trees_count >= max ) {
                        break;
                    }
                    tree_index++;
                }
                else if( String_i_start_with(reader->line, "END", true)){
                    break;
                }
            }
        }
    }
    
    free_FileReader(reader);
    free_StringBuffer(buffer);
    
    //fprintf(stderr, "count %d\n", *trees_count);
    
    if ( names_count != 0 ) {
        for (int i = 0; i < *trees_count; i++) {
            trees[i] = insertNexusNames(trees[i], names);
        }
        for (int i = 0; i < names_count; i++) {
            free(names[i]); names[i] = NULL;
        }
        free(names); names = NULL;
    }
    
    
    
    if ( positions != NULL){
        for ( int i = 0; i < nRange; i++) {
            free(positions[i]);
        }
        free(positions);
    }
    
    if ( *trees_count == 0 ) {
        free(trees);
        trees = NULL;
    }
    return trees;
}


// Read newick file but read only the first tree
char *readNewickTree( const char *infile ){
    char *nexus = NULL;
    
    FileReader *reader = new_FileReader(infile, 1000);
    
    int count = 0;
    while ( reader->read_line(reader) ) {
        if ( count == 0 && String_i_start_with(reader->line, "(", true) ) {
            char *ptr = &reader->line[ strlen(reader->line)-1 ];
            while( *ptr != ';' && *ptr != *reader->line ){
                ptr--;
            }
            if ( *ptr == *reader->line ) {
                error("missing ; in tree file\n");
            }
            *ptr = '\0';
            
            ptr = reader->line;
            while ( *ptr != '(') {
                ptr++;
            }
            nexus = String_clone( ptr );
            count++;
        }
        else if ( count != 0 && String_i_start_with(reader->line, "(", true) ) {
            fprintf(stderr, "Only the first tree is read, the rest is ignored\n");
            break;
        }
    }
    
    free_FileReader(reader);
    return nexus;
}

char **readNewickTrees( const char *infile, int *trees_count ){
    int trees_capcacity = 10;
    char **trees = (char**)calloc( trees_capcacity, sizeof(char*) );
    
    FileReader *reader = new_FileReader(infile, 1000);
    
    char *line = NULL;
    
    int count = 0;
    while ( reader->read_line(reader) ) {
        line = reader->line;
        line = String_trim(line);
        
        if( strlen(line) == 0 ) continue;
        

        char *ptr = &reader->line[ strlen(reader->line)-1 ];
        while( *ptr != ';' && *ptr != *reader->line ){
            ptr--;
        }
        if ( *ptr == *reader->line ) {
            error("missing ; in tree file\n");
        }
        *ptr = '\0';
        
        ptr = reader->line;
        while ( *ptr != '(') {
            ptr++;
        }
        
        if( count == trees_capcacity){
            trees_capcacity += 10;
            trees = realloc(trees, trees_capcacity*sizeof(char*));
        }
        
        trees[count] = String_clone( ptr );
        count++;
        
    }
    assert(count>0);
    if( trees_capcacity != count ){
        trees = realloc(trees, count*sizeof(char*));
    }
    free_FileReader(reader);
    *trees_count = count;
    return trees;
}


#pragma mark -
// MARK: RESCUE

char * readGreedyLogTree( const char *infile, const double alpha, double *lnl, int *failure ){
	
	const char *TREE_TAG = "TREE";
	const char *NEXUS_TREE_BLOCK = "BEGIN TREES";
	const char *NEXUS_TREE_TRANSLATE_BLOCK = "TRANSLATE";
	
	
	unsigned names_count = 0;
	char **names = NULL;
	
	FileReader *reader = new_FileReader(infile, 1000);
	StringBuffer *buffer = new_StringBuffer(200);
	
	int nLocalClocks_previous = 1;
	int tree_count_previous = 0;
	double lk_max_previous = -10000000000;
	StringBuffer *best_tree_previous = new_StringBuffer(500);
	
	int nLocalClocks_current = 1;
	int tree_count_current = 0;
	double lk_max_current = -10000000000;
	StringBuffer *best_tree_current = new_StringBuffer(500);
	
	bool done = false;
	
	while ( reader->read_line(reader) ) {
		
		// PARSE TREE BLOCK
		if ( String_i_start_with(reader->line, NEXUS_TREE_BLOCK, true) ) {
			while ( reader->read_line(reader) ) {
				// PARSE TRANSLATE BLOCK
				if ( String_i_start_with(reader->line, NEXUS_TREE_TRANSLATE_BLOCK, true) ) {
					names = parseNexusTranslateBlock(reader, buffer, &names_count);
				}
				
				
				
				if ( String_i_start_with(reader->line, TREE_TAG, true) ) {
					char *ptr = reader->line;
					
					// get likelihood
					char *pch = strstr (ptr,"lk=");
					pch += 3;
					StringBuffer_empty(buffer);
					while ( *pch != ']' ) {
						StringBuffer_append_char(buffer, *pch);
						pch++;
					}
					double lk = atof(buffer->c);
					
					// get number of local clocks
					pch = strstr (ptr,"local=1");
					int nLocalClocks = 0;
					while ( pch != NULL ) {
						nLocalClocks++;
						pch += 5;
						pch = strstr (pch,"local=1");
					}
					
					
					// check if the newcik tree is not truncated (in case of a crash)
					int pos = strrchr(ptr,';') - ptr;
					if ( pos < 0 ) {
						goto END;
					}
					
					if ( nLocalClocks != nLocalClocks_current ) {
						if ( nLocalClocks_current == 1 ) {
							StringBuffer_empty(best_tree_previous);
							StringBuffer_append_StringBuffer(best_tree_previous, best_tree_current);
							lk_max_previous = lk_max_current;
							tree_count_previous = tree_count_current;
						}
						else {
							double p = LRT(lk_max_previous, lk_max_current, nLocalClocks_previous, nLocalClocks_current );
							fprintf(stderr, "LnL1: %f LnL0: %f p: %f df: %d - %d ", lk_max_current, lk_max_previous, p, nLocalClocks_current, nLocalClocks_previous);
							if ( p < alpha) {
								StringBuffer_empty(best_tree_previous);
								StringBuffer_append_StringBuffer(best_tree_previous, best_tree_current);
								lk_max_previous = lk_max_current;
								nLocalClocks_previous = nLocalClocks_current;
								tree_count_previous = tree_count_current;
								fprintf(stderr, "*\n");
							}
							else {
								fprintf(stderr, "\n");
								done = true;
								goto END;
							}
							
						}
						nLocalClocks_current = nLocalClocks;
						tree_count_current = 0;
						
					}
					
					if( lk > lk_max_current ){
						lk_max_current = lk;
						bool comment = false;
						while ( *ptr != '(' || comment ) {						
							if ( *ptr == '[' ) {
								comment = true;
							}
							else if( comment && *ptr == ']' ) {
								comment = false;
							}
							ptr++;
						}
						
						int pos = strrchr(ptr,';') - ptr;
						if( pos >= 0) ptr[pos] = '\0';
						else break;
						
						StringBuffer_set_string(best_tree_current, ptr);
					}
					tree_count_current++;
					
					
				}
				else if( String_i_start_with(reader->line, "END", true)){
					break;
				}
			}
		}
	}
END:
	free_FileReader(reader);
	
	char *tree = NULL;
	//int nClocks = 0;
	
	if ( tree_count_previous-1 != tree_count_current  || done ) {
		tree = best_tree_previous->c;
		//nClocks = nLocalClocks_previous;
		*lnl = lk_max_previous;
		free(best_tree_previous);
		free_StringBuffer(best_tree_current);
		fprintf(stderr, "tree count previous %d current %d\n", tree_count_previous, tree_count_current);
		*failure = 0;
	}
	// we have reached the end of the file and no better model was found, try the last one
	else {
		double p = LRT(lk_max_previous, lk_max_current, nLocalClocks_previous, nLocalClocks_current );
		fprintf(stderr, "LnL1: %f LnL0: %f p: %f df: %d - %d ", lk_max_current, lk_max_previous, p, nLocalClocks_current, nLocalClocks_previous);
		if ( p < alpha) {
			tree = best_tree_current->c;
			//nClocks = nLocalClocks_current;
			*lnl = lk_max_current;
			free(best_tree_current);
			free_StringBuffer(best_tree_previous);
			*failure = 1;
			fprintf(stderr, "*\n");
		}
		else {
			tree = best_tree_previous->c;
			//nClocks = nLocalClocks_previous;
			*lnl = lk_max_previous;
			free(best_tree_previous);
			free_StringBuffer(best_tree_current);
			*failure = 0;
			fprintf(stderr, "\n");
		}
		
	}
	
	if ( names_count != 0 ) {
		tree = insertNexusNames( tree, names );
		
		for (int i = 0; i < names_count; i++) {
			free(names[i]);
		}
		free(names);
	}
	
	//fprintf(stderr, "Number of local clocks: %d LnL: %f\n", nClocks, *lnl );
	
	free_StringBuffer(buffer);
	
	return tree;
}
