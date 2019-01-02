/*
 *  tree.h
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

#ifndef _TREE_H_
#define _TREE_H_

#include "utils.h"
#include "parameters.h"
#include "node.h"
#include "mjson.h"
#include "nj.h"
#include "upgma.h"

typedef enum  treeorder {PREORDER, POSTORDER} treeorder;

struct _Tree;
typedef struct _Tree Tree;


#pragma mark -
#pragma mark Tree




Model * new_TreeModel( const char* name, Tree *tree );

Model* new_TreeModel_from_json(json_node* node, Hashtable* hash);

Tree * new_Tree( const char *nexus, bool containBL );

Tree * new_Tree2( Node *root, bool containBL );

void Tree_init_heights ( Tree *atree );

void free_Tree( Tree *t);

Tree * clone_Tree( const Tree *tree );

Tree * clone_SubTree( const Tree *tree, Node *node );

int Tree_node_count( const Tree *tree );

int Tree_tip_count( const Tree *tree );

Node * Tree_root( Tree *tree );

void Tree_set_root( Tree *tree, Node *root );

int Tree_id( const Tree *tree );

void Tree_set_id( Tree *tree, int id );

void Tree_set_dated( Tree *tree, bool dated );

bool Tree_dated( const Tree *tree );

void create_node_list( Tree *atree, treeorder order );

Node ** get_tips( Tree *atree, treeorder order );

Node ** Tree_get_nodes( Tree *atree, treeorder order );

Node * Tree_get_node(Tree *atree, treeorder order, int index);

Node * Tree_get_node_by_name(Tree *tree, const char *name );

Node ** Tree_nodes( Tree *tree );

Node * Tree_node( Tree *tree, int index );


void postorder_generic( Tree *tree, void *data, void (*action)( Node *node, void *data ));
void preorder_generic( Tree *tree, void *data, void (*action)( Node *node, void *data ));

void Tree_set_topology_changed( Tree *tree );

bool Tree_topology_changed( const Tree *tree );

void Tree_update_topology( Tree *tree );

void rename_tree( Tree *tree, treeorder order );

void Tree_init_depth( Tree *tree );

void Tree_set_rooted(Tree* tree, bool rooted);

bool Tree_rooted(Tree* tree);

#pragma mark -
#pragma mark Height
void calculate_heights( Tree *tree );

bool parse_dates( Tree *tree );

bool Tree_is_calibrated( Tree *tree );

bool parse_calibrations( Tree *tree );

void Tree_init_heights_heterochronous( Tree *tree, const double rate, bool convert );

void Tree_init_heights_homochronous( Tree *tree, const double rate );
									
void summarize_calibrations( Tree *tree );

double mean_rate_calibrations( Tree *tree );

void Tree_constrain_height( Node *node );

void Tree_constraint_heights( Tree *tree );

void Tree_constrain_neighbor_heights( Node *node );

void Tree_check_constraint_heights( Tree *tree );

void Tree_scale_heights( Node *node, const double scaler );


void Tree_heights_to_vector( Tree *tree, double *heights );
void Tree_vector_to_heights( const double *heights, Tree *tree );

void Tree_restore_branch_length( Tree *tree, const Parameters *bl );

double Tree_scale_total_length( Tree *tree, double total );

void Tree_branch_length_to_vector( Tree *tree, double *distances );

void Tree_vector_to_branch_length( Tree *tree, const double *distances );

void Tree_copy_heights( Tree *src, Tree *dist );

// Print

void print_compare_tree(  Tree *t1,  Tree *t2);

void compare_tree( const Tree *t1, const Tree *t2 );

void Tree_print_parameters( Tree *tree );

double Tree_distance_to_root( const Node *node );

double Tree_distance_between_nodes( const Node *node, const Node *ancestor );



void Tree_newick_distance( StringBuffer *buffer, Node *n );

void Tree_to_newick( StringBuffer *buffer, Node *n );

void Tree_heights_to_times( Node *node, double *times );

void Tree_set_prop_height( Tree *tree, const double value, const int index );


//void Tree_scale2( Tree *tree, double scale );

//void Tree_scale( Tree *tree, double scale );

void Tree_scale_distance( Tree *tree, double scale );

//void Tree_shift_up( Tree *tree, const double shift );

void Tree_expand( Tree *tree, const double shift );

//void Tree_slide_node_up( Tree *tree, Node *node, double amount );

//void Tree_slide_node_down( Tree *tree, Node *node, double amount );


void Tree_reroot( Tree *tree, Node *node );

void Tree_rearrange( Tree *tree, Node *node1, Node *node2 );

void Tree_swap_parents_by_name( Tree *tree, char *name1, char *name2 );

void Tree_post_to_preorder( Tree *tree, double *v );

void Tree_copy_distances( const Tree* src, Tree *dst );

void Tree_copy_height_constraints( Tree *src, Tree *dist );


double Tree_length( const Tree *tree );


char * Tree_to_string_nexus_with_annotation( Tree *tree );

void Tree_StringBuffer_nexus_with_annotation( StringBuffer *treebuffer, Tree *tree );

bool Tree_is_time_mode(Tree* tree);

#endif

