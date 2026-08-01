// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_node_h
#define PhyC_node_h

#include "parameters.h"
#include "hashtable.h"

#define BL_DEFAULT 0.1

#define BL_MIN 1.E-8
#define BL_MAX 10


#define POSTFIX_DISTANCE "distance"
#define POSTFIX_HEIGHT   "height"
#define CALIBRATION_TAG  "cal_height"

// Sentinel for Node.branch_index: the node owns no entry in the distance vector.
#define NODE_NO_BRANCH ((size_t)-1)

typedef struct Node{
	int id;
    char *name;
    struct Node *parent;
    struct Node *right;
    struct Node *left;
    int depth;
	
	double time; // contains age of node as given by the user
	
	int preorder_idx;
	int postorder_idx;
	
	Parameter *distance;
	// Index of this node's branch in the shared distance vector. It is NOT the
	// node id: an unrooted tree stored as a binary tree has 2n-2 branch-owning
	// nodes but only 2n-3 free branch lengths, so the two children of the root
	// share one entry. NODE_NO_BRANCH for the root and for the child of the root
	// whose branch is pinned to zero. Maintained by Tree_update_branch_indices.
	size_t branch_index;
	Parameter *height; // constraint lower bound maximum of the 2 kids

	int class_id;
	
	char *info;
	double bl;// temporary value for branch length when parsing nexus
    Hashtable *annotation;
	bool poly;
} Node;

Node * new_EmptyNode();

Node * new_Node( Node *parent, const char *nodename, const int counter );

void free_Node( Node *node );

Node * clone_Node( const Node *node);

void Node_set_height( Node *node, const double value );
void Node_set_height_quietly( Node *node, const double value );

double Node_height( const Node *node );

double Node_time_elapsed( Node *node );

void Node_set_distance( Node *node, const double value );

double Node_distance( const Node *node );

bool Node_has_distance( const Node *node );

void Node_set_parent( Node *node, Node *parent );

Node * Node_parent( Node *node );

bool Node_isroot( const Node *node );

bool Node_isleaf( const Node *node );

Node * Node_left( Node *node );

Node * Node_right( Node *node );

Node * Node_sibling( Node *node );

bool Node_isancestor(  Node  * const node, const Node *putative_ancestor );

bool Node_isrelated( Node *aNode1, Node *aNode2 );

int Node_graph_distance( Node *node1, Node *node2 );

void Node_swap_parents( Node *node1, Node *node2 );

void Node_rotate( Node *node );

int Node_tip_count( Node *node );

double Node_time( Node *node );

void Node_set_time( Node *node, const double time );

int Node_id( const Node *node );

void Node_set_id( Node *node, const int id );

int Node_class_id( const Node *node );

void Node_set_class_id( Node *node, const int class_id );

int Node_depth( Node *node );

void Node_set_depth( Node *node, const int depth );

char * Node_name( Node *node );

void Node_set_name( Node *node, const char *name );

void Node_rename( Node *n, const int newPos );

void Node_set_annotation( Node *node, const char *key, const char *value );

void Node_empty_annotation( Node *node );

char * Node_annotation( Node *node, const char *key );

double Node_get_double_from_info( const Node *node, const char *str );

int Node_get_int_from_info( const Node *node, const char *str );

char * Node_get_string_from_info( const Node *node, const char *str );

void Node_removeChild( Node *parent, Node *node);

bool Node_addChild( Node *parent, Node *node);

void Node_add_listener( Node *node, Model *model );

#endif
