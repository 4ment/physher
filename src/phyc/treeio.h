/*
 *  treeio.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/1/10.
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

#ifndef _TREEIO_H_
#define _TREEIO_H_
#include <stdio.h>

#include "tree.h"
#include "filereader.h"

struct _TreeFileIterator;
typedef struct _TreeFileIterator TreeFileIterator;

struct _TreeFileIterator {
    FileReader *reader;
    Hashtable *translation;
	bool more;
};



#pragma mark Iterator

TreeFileIterator *new_TreeFileIterator( const char *filename );

void free_TreeFileIterator(TreeFileIterator *iter );

char * TreeFileIterator_next_tree( TreeFileIterator *iter );

#pragma mark -

char *readTree( const char *infile );

char **readTrees( const char *infile, int *trees_count );

char **readTreesRange( const char *infile, int *trees_count, int start, int end );

// Nexus
char ** readNexusTreesRange( const char *infile, int *trees_count, int start, int end );

char ** readNexusTrees( const char *infile, int *trees_count, char *range, int max );

char * readNexusFirstTree( const char *infile );

char * readNexusTreeAt( const char *infile, const int index );

char ** readNexusNFirstTree( const char *infile, int *count, int max );

char ** readNexusRangeTree( const char *infile, int *count, char *range );


char ** parseNexusTranslateBlock( FileReader *reader, StringBuffer *buffer, unsigned *count );
char * insertNexusNames( char *tree, char ** const names );
void insertNexusNamesStringBuffer( StringBuffer *tree, char ** const names );

// Newick tree
char *readNewickTree( const char *infile );

char **readNewickTrees( const char *infile, int *trees_count );

#pragma mark -
#pragma mark Print

// Newick
void Tree_print_newick( FILE *pf, Tree *tree, bool internal );

void Tree_print_newick_subtree( FILE *pf, bool time, const Node *n, bool internal );

void Tree_print_height_newick( FILE *pf, Tree *tree, bool internal );

void print_tree_extended( FILE *pf, const Node *n, char **info );

void print_tree_extended_nexus_double( FILE *pf, Tree *tree, double *info );


// Nexus

//void print_nexus_header_figtree( FILE *pf, Tree *t, const char **orderedNames );
//void print_nexus_header_figtree_Taxa( FILE *pf, Tree *tree, const char **orderedNames );
//void print_nexus_header_figtree_BeginTrees( FILE *pf, Tree *tree, const char **orderedNames );
//void Tree_print_with_info( FILE *pf, Tree *tree, const Node *n, const char **names );


// These functions index each tip with its preorder index
// Should not be used with trees with different topologies
void Tree_print_nexus_header_figtree( FILE *pf, Tree *tree );
void Tree_print_nexus_header_figtree_Taxa( FILE *pf, Tree *tree );
void Tree_print_nexus_taxa_block( FILE *pf, Tree *tree );
void Tree_print_nexus_header_figtree_BeginTrees( FILE *pf, Tree *tree );
void Tree_print_nexus_with_annotation2( FILE *pf, Tree *tree, bool time );
void Tree_print_nexus_with_annotation( FILE *pf, Tree *tree );
void Tree_print_nexus( FILE *pf, Tree *tree );



// RESCUE
char * readGreedyLogTree( const char *infile, const double alpha, double *lnl, int *failure );

#endif
