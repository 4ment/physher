/*
 *  tree.c
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

#include "tree.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>

#include "mstring.h"
#include "matrix.h"
#include "mathconstant.h"
#include "node.h"

#include "distancematrix.h"

struct _Tree{
	int id;
	Node *root;
	int nTips;
	int nNodes;
	bool rooted;
	Node **postorder;
	Node **preorder;
    Node **nodes;
	bool topology_changed;
	bool dated;
	bool time_mode;
	Parameters* distances;
};

static void tree_init_distance_parameters( Tree *t, Parameters *ps );
static void tree_init_height_parameters( Tree *t, Parameters *ps );
static void _Tree_count_nodes( Node *node, int *tips, int *nodes);


void check_tree( Tree *tree ){
	Node ** nodes = Tree_get_nodes(tree, POSTORDER);
	for ( int i = 0; i < tree->nNodes; i++) {
		if(tree->root != nodes[i]) assert( nodes[i]->parent);
	}
}

#pragma mark -
#pragma mark Tree

// called just after a ) or a leaf name
// can start  with bootstrap value, node name, :, [
// return i on ')' or ',' or ;
void _parse_description( Node *node, const char *newick, int *i,  bool containBL ){
	
	while ( isspace( newick[*i] ) ) ++(*i);
	
	if( newick[*i] == ',' || newick[*i] == ')' || newick[*i] == ';' ) return;

	StringBuffer *buffer = new_StringBuffer(100);
	
	double rate   = -INFINITY;
	double height = -INFINITY;
	double d      = -INFINITY;
	
	if( newick[*i] == '[' ){
		do {
			StringBuffer_append_char(buffer, newick[*i]);
			++(*i);			
		} while ( newick[*i] != ']');
		StringBuffer_append_char(buffer, newick[*i]);
		
		node->info = StringBuffer_assign(buffer, node->info);
		++(*i);
	}
		
	// bootstrap or internal node name are ignored for now
	while ( newick[*i] != ':' && newick[*i] != ',' && newick[*i] != ')' && newick[*i] != ';' ) {
		++(*i);	
	}
	
	
	if( newick[*i] == ':' ){
		StringBuffer_empty(buffer);
		(*i)++;
		// This will remove the first [&...]
		if( newick[*i] == '[' ){
			do {
				StringBuffer_append_char(buffer, newick[*i]);
				++(*i);
				
			} while ( newick[*i] != ']');
			StringBuffer_append_char(buffer, newick[*i]);
			
			node->info = StringBuffer_assign(buffer, node->info);
			++(*i);
		}
		
		StringBuffer_empty(buffer);
		
		while ( newick[*i] != ',' && newick[*i] != ')' ){
			StringBuffer_append_char(buffer, newick[*i] );
			(*i)++;
		}
		if(buffer->length == 0) fprintf(stderr, "info: %s -%s-\n", Node_name(node), buffer->c);
		d = atof(buffer->c);
	}
	
	
	if ( node->info != NULL ) {
		if ( String_contains_str(node->info, "rate=")) {
			char *pch = node->info + String_index_of_str(node->info, "rate=") + strlen("rate=");
			StringBuffer_empty(buffer);
			while ( *pch != ',' && *pch != ']' ){
				StringBuffer_append_char(buffer, *pch );
				pch++;
			}
			
			rate = atof(buffer->c);
			
			if ( String_contains_str(node->info, "height=")) {
				pch = node->info + String_index_of_str(node->info, "height=") + strlen("height=");
				StringBuffer_empty(buffer);
				while ( *pch != ',' && *pch != ']' ){
					StringBuffer_append_char(buffer, *pch );
					pch++;
				}
				
				height = atof(buffer->c);
			}
		}
	}
	
	if ( containBL ) {
		Node_set_distance(node, dmax(BL_MIN, d));
		if( height != -INFINITY ){
			Node_set_height(node, height);
		}
		if( rate != -INFINITY ){// that's wrong
			Node_set_height(node, Node_distance(node) / rate );
		}
	}
	else {
		if( d != -INFINITY ){
            Node_set_height(node, d);
            if( rate != -INFINITY ){
                Node_set_distance(node, d*rate);
                //printf("%f %f %f\n",d,rate, (d*rate));
            }
            else {
                Node_set_distance(node, d);
            }
        }
	}
	
	//fprintf(stderr, "%s %s\n", node->name, (node->info == NULL ? "no" : node->info));
	free_StringBuffer(buffer);
}

Tree * new_Tree( const char *nexus, bool containBL ){
	Tree *atree = (Tree *)malloc(sizeof(Tree));
	assert(atree);
	atree->root = NULL;
	atree->postorder = NULL;
	atree->preorder = NULL;
	atree->nodes = NULL;
	atree->id = 0;
	atree->nTips = 0;
	atree->nNodes = 0;
	atree->rooted = false;
	atree->dated = false;
	atree->time_mode = false;

	//printf("%s",nexus);
	Node *current = NULL;
	int i = 0;
	unsigned l = strlen(nexus);
	
	StringBuffer *buffer = new_StringBuffer(50);
	
	int count = 0;
	for ( i = 0; i < l; i++) {
		if ( nexus[i] == '(' ) {
			count++;
		}
		else if ( nexus[i] == ')' ) {
			count--;
		}
	}
	
	if ( count != 0 ){
		fprintf(stderr, "The newick tree is malformed: Number opening parenthesis != number closing parenthesis\n");
		exit(1);
	}
    
    //assert(nexus[0] == '(');
    current = atree->root = new_Node(NULL, NULL, atree->nNodes);
    atree->nNodes++;
    
	for ( i = 1; i < l; i++) {
		// Internal node
		if( nexus[i] == '(' ){
			Node *n = new_Node(current, NULL, atree->nNodes);
            
            bool success = Node_addChild(current, n);
            // polytomy
            if( !success ) {
                Node *temp = new_Node(current, NULL, atree->nNodes-1);
                temp->poly = true;
                
                Node_set_distance(temp, BL_MIN);
                
                Node *r = current->right;
                current->right = temp;
                temp->left = r;
                temp->right = n;
                n->parent = temp;
                r->parent = temp;
                Node_rename(n, atree->nNodes);
                atree->nNodes++;
            }
			
			
			atree->nNodes++;
			current = n;
		}
		else if( nexus[i] == ',' ){
			//fprintf(stderr, "%c\n", nexus[i]);
		}
		else if( nexus[i] == ':' ){
			fprintf(stderr, "%s\n", &nexus[i] );
		}
		else if( nexus[i] == ')' ){
			// we have reached the end and no info after the last )
			if( atree->root == current && i == l-1 ){
                
				
			}
			// we have reached the last ) but there is some info for the root
			else if(atree->root == current ){
                //fprintf(stderr, "-- %s\n", &nexus[i] );
				i++;
				StringBuffer_empty(buffer);
				while ( i != l ){
					StringBuffer_append_char(buffer, nexus[i]);
					i++;
				}
				
				bool isbootstrap = isFloat(buffer->c);
				//fprintf(stderr, "bootstartp %lu =%s=, -%s-\n",buffer->length, buffer->c,&nexus[i]);
				
				if ( isbootstrap ) {
					StringBuffer_prepend_string(buffer, "[&default=" );
					StringBuffer_append_string(buffer, "]");
					if ( current->info == NULL ) {
						current->info = String_clone( buffer->c );
					}
					else {
						current->info = String_append_string(current->info, buffer->c);
					}
				}
                else {
                    //_parse_description(current, nexus, &i, containBL);
                    //printf(" ---- %s\n", buffer->c);
                    current->info = String_clone( buffer->c );
                }
				
				// ignore for now
				if ( buffer->c[0] == ':' ) {
					//i++;
					//getInfo(current, nexus, &i, containBL );
				}
				
			}
			// internal node
			else {
				i++;
				_parse_description(current, nexus, &i, containBL);
				i--;
				
				current = Node_parent(current);
				if ( current->poly == true ) {
					current = Node_parent(current);
				}
			}
			
		}
		// Node ID: tip or internal
		// does not parse bootstrap or other values at each node
		else{
			StringBuffer_empty(buffer);
			// read name
			while ( nexus[i] != ':' && nexus[i] != ',' && nexus[i] != ')' && nexus[i] != '[' ){
				StringBuffer_append_char(buffer, nexus[i]);
				i++;
			}
			assert( buffer->length != 0 );
			char *ptr = buffer->c;
			if ( (buffer->c[0] == '\'' && buffer->c[buffer->length-1] == '\'')
				|| (buffer->c[0] == '"' && buffer->c[buffer->length-1] == '"') ) {
				buffer->c[buffer->length-1] = '\0';
				ptr++;
			}
			
			Node *n = new_Node(current, ptr, atree->nNodes);
			
			_parse_description(n, nexus, &i, containBL);
			--i;
			
			assert(current);
			if(current->left == NULL ){
				current->left = n;
				current = Node_parent(n);
			}
			else if(current->right == NULL ){
				current->right = n;
				current = Node_parent(n);
			}
			else {
				Node *temp = new_Node(current, NULL, atree->nNodes-1);
				temp->poly = true;
				Node_set_distance(temp, 0.0);
				
				//fprintf(stderr, "%s %s n:%s\n", current->name, temp->name, n->name);
				
				Node *r = current->right;
				current->right = temp;
				temp->left  = r;
				temp->right = n;
				n->parent = temp;
				r->parent = temp;
				
				atree->nNodes++;
				//fprintf(stderr, "poly2\n");
				current = Node_parent(temp);
			}
			
			atree->nTips++;
			atree->nNodes++;
		}
	}
	
	//assert( atree->root == current );
	if( atree->root != current ){
		fprintf(stderr, "%s %s\ncurrent different from root\n", atree->root->name, current->name );
		exit(1);
	}
	
	Node_set_distance( current, -1 );
	//Parameter_set_bounds( current->distance, -1, -1 );
	//Parameter_set_fixed( current->distance, true );
	
	Tree_update_topology(atree);
    
    Node **nodes = Tree_get_nodes(atree, POSTORDER);
    atree->nodes = (Node**)malloc( atree->nNodes* sizeof(Node*) );
    assert(atree->nodes);
    
    int nTips = 0;
    int nInternals = 0;
    for ( int i = 0; i < Tree_node_count(atree); i++ ) {
        nodes[i]->id = i;
		nodes[i]->distance->id = i;
		nodes[i]->height->id = i;
        if( Node_isleaf(nodes[i]) ){
            nodes[i]->class_id = nTips++;
        }
        else {
            nodes[i]->class_id = nInternals++;
        }
        atree->nodes[ nodes[i]->id ] = nodes[i];
    }
	
	check_tree(atree);
	
	// is going to fail for homochronous trees
	if ( !containBL ) {
		Tree_init_heights(atree);
	}
	
	atree->distances = new_Parameters(Tree_node_count(atree)-2);
	Node* root = Tree_root(atree);
	for (int i = 0; i < Tree_node_count(atree); i++) {
		Node* node = Tree_node(atree, i);
		if ( node != root && root->right != node) {
			Parameters_add(atree->distances, node->distance);
		}
	}
	
	free_StringBuffer(buffer);
	
	return atree;
}

Tree * new_Tree2( Node *root, bool containBL ){
	Tree *atree = (Tree *)malloc(sizeof(Tree));
	assert(atree);
	atree->root = root;
	atree->postorder = NULL;
	atree->preorder = NULL;
    atree->nodes = NULL;
	atree->id = 0;
	atree->nTips = 0;
	atree->nNodes = 0;
	atree->rooted = false;
	atree->dated = false;
	atree->time_mode = false;
	
	_Tree_count_nodes(root, &atree->nTips, &atree->nNodes);
	
	Node_set_distance( root, -1 );
	//Parameter_set_bounds( root->distance, -1, -1 );
	//Parameter_set_fixed( root->distance, true );
	
	Tree_update_topology(atree);
    
    Node **nodes = Tree_get_nodes(atree, POSTORDER);
    atree->nodes = (Node**)malloc( atree->nNodes* sizeof(Node*) );
    assert(atree->nodes);
    int nTips = 0;
    int nInternals = 0;
    for ( int i = 0; i < Tree_node_count(atree); i++ ) {
		nodes[i]->id = i;
		nodes[i]->distance->id = i;
		nodes[i]->height->id = i;
        if( Node_isleaf(nodes[i]) ){
            nodes[i]->class_id = nTips++;
        }
        else {
            nodes[i]->class_id = nInternals++;
        }
        atree->nodes[ nodes[i]->id ] = nodes[i];
    }
	
	check_tree(atree);
	
	// is going to fail for homochronous trees
	if ( !containBL ) {
		Tree_init_heights(atree);
	}
	
	atree->distances = new_Parameters(Tree_node_count(atree)-2);
	for (int i = 0; i < Tree_node_count(atree); i++) {
		Node* node = Tree_node(atree, i);
		if ( node != root && root->right != node) {
			Parameters_add(atree->distances, node->distance);
		}
	}
	return atree;
}



// This function is called when the tree branches contain time in the newick file
void Tree_init_heights ( Tree *atree ) {

	parse_dates(atree);

	double *heights = dvector( Tree_node_count(atree) );
	bool *set = bvector( Tree_node_count(atree) );
	Tree_heights_to_vector(atree, heights);
	
	// heights in the nexus file have to be transformed
	Node **nodes = Tree_get_nodes(atree, POSTORDER);
	
	// if the sequences are homochronous, time is set to 0
	bool need_conv = true;
	double youngest = 0;
	
	
	for ( int i = 0; i < Tree_node_count(atree); i++) {	
		if ( Node_isleaf(nodes[i]) ){
			if( nodes[i]->time == 0.0 ) need_conv = false;
			youngest = dmax(youngest, Node_time(nodes[i]));
		}				
	}
	
	
	//if(need_conv)fprintf(stderr, "need conversion\n");
	for ( int i = 0; i < Tree_node_count(atree); i++) {
		if( Node_isleaf(nodes[i]) ){
			nodes[i]->height->value = (need_conv ? youngest - nodes[i]->time : nodes[i]->time );
			assert(nodes[i]->height->value >=0 );
		}
	}
	
	StringBuffer *buffer = new_StringBuffer(10);
	for ( int i = 0; i < Tree_node_count(atree); i++) {
		//fprintf(stderr, "init %s %s\n", Node_name(Node_parent(nodes[i])), Node_name(nodes[i]) );
		
		if( !Node_isroot(nodes[i]) && !set[Node_parent(nodes[i])->postorder_idx] ){
			Node_set_height( Node_parent(nodes[i]), Node_height(nodes[i]) + heights[i] );
            //printf("nonode %s %f %f %E\n", nodes[i]->name, Node_parent(nodes[i])->height->value, nodes[i]->height->value, Node_parent(nodes[i])->height->value- nodes[i]->height->value);
			assert(Node_parent(nodes[i])->height->value >=0 );
			//assert(Node_parent(nodes[i])->height->value > nodes[i]->height->value);
			set[Node_parent(nodes[i])->postorder_idx] = true;
		}
		
		// Grab calibration points
		if ( nodes[i]->info != NULL ) {
            //printf("%s %s\n", nodes[i]->name, nodes[i]->info);
			int pos = String_index_of_str(nodes[i]->info, CALIBRATION_TAG);
			if ( pos != -1 ) {
				char *ptr = nodes[i]->info+pos;
				while ( *ptr != '{' ) {ptr++;}
				ptr++;
				StringBuffer_empty(buffer);
				
				while ( *ptr != ',' ) {
					StringBuffer_append_char(buffer, *ptr);
					ptr++;
				}
				
				if ( buffer->length > 0 ) {
					double lower = atof(buffer->c);
					Parameter_set_lower(nodes[i]->height, lower );
					Constraint_set_flower(nodes[i]->height->cnstr, lower);
					Constraint_set_lower_fixed(nodes[i]->height->cnstr, true);
				}
				
				StringBuffer_empty(buffer);
				
				ptr++;
				
				while ( *ptr != '}' ) {
					StringBuffer_append_char(buffer, *ptr);
					ptr++;
				}
				if ( buffer->length > 0 ) {
					double upper = atof(buffer->c);
					Parameter_set_upper(nodes[i]->height, upper);
					Constraint_set_fupper(nodes[i]->height->cnstr, upper);
					Constraint_set_upper_fixed(nodes[i]->height->cnstr, true);
				}
				
				if ( Constraint_lower_fixed(nodes[i]->height->cnstr) && Constraint_upper_fixed(nodes[i]->height->cnstr) && Parameter_upper(nodes[i]->height) == Parameter_lower(nodes[i]->height) ) {
					Parameter_set_estimate(nodes[i]->height, false);
				}
				
			}
		}
		
		
	}
	
	//infer_distance_from_rate_height(bm);
	
	Tree_constraint_heights(atree);
	free(heights);
	free(set);
	free_StringBuffer(buffer);
}

void _tree_handle_change( Model *self, Model *model, int index ){
	Tree *tree = (Tree*)self->obj;
	if ( tree->time_mode ) {
		Node *node = Tree_node(tree, index);
		// constrain lower bound of parent
		if( !Node_isroot(node) ) {
			Tree_constrain_height( Node_parent(node) );
		}
		
		// constrain upper bound of children
		if( !Node_isleaf(node) ){
			Tree_constrain_height( Node_left(node) );
			Tree_constrain_height( Node_right(node) );
		}
		
		self->listeners->fire( self->listeners, self, Node_id(Node_left(node)) );
		self->listeners->fire( self->listeners, self, Node_id(Node_right(node)) );
	}
	self->listeners->fire( self->listeners, self, index );
}

static void _tree_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free tree model %s\n", self->name);
		Tree *tree = (Tree*)self->obj;
		free_Tree(tree);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _tree_model_clone( Model *self, Hashtable *hash ){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	Tree *tree = (Tree*)self->obj;
	Tree *clonetree = clone_Tree(tree);
	for (int i = 0; i < Tree_node_count(clonetree); i++) {
		Node* node = Tree_node(clonetree, i);
		if(node->distance != NULL){
			char* name = Parameter_name(node->distance);
			Hashtable_add(hash, name, node->distance);
		}
		if(node->height != NULL){
			char* name = Parameter_name(node->height);
			Hashtable_add(hash, name, node->height);
		}
	}
	Model* clone = new_TreeModel(self->name, clonetree);
	Hashtable_add(hash, clone->name, clone);
	return clone;
}

static void _tree_model_get_free_parameters(Model* model, Parameters* parameters){
	Tree* mtree = (Tree*)model->obj;
	for ( int i = 0; i < Tree_node_count(mtree); i++ ) {
		Node* node = Tree_node(mtree, i);
		if(Parameter_estimate(node->distance)){
			Parameters_add(parameters, node->distance);
		}
		if(Parameter_estimate(node->height)){
			Parameters_add(parameters, node->height);
		}
	}
}

// TreeModel listen to the height and distance parameters
Model * new_TreeModel( const char* name, Tree *tree ){
	Model *model = new_Model("tree", name, tree);
	StringBuffer* buffer = new_StringBuffer(10);
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		StringBuffer_set_string(buffer, name);
		StringBuffer_append_strings(buffer, 2, ".", tree->nodes[i]->distance->name);
		Parameter_set_name(tree->nodes[i]->distance, buffer->c);
		StringBuffer_set_string(buffer, name);
		StringBuffer_append_strings(buffer, 2, ".", tree->nodes[i]->height->name);
		Parameter_set_name(tree->nodes[i]->height, buffer->c);
		tree->nodes[i]->distance->listeners->add(tree->nodes[i]->distance->listeners, model);
		tree->nodes[i]->height->listeners->add(tree->nodes[i]->height->listeners, model);
	}
	free_StringBuffer(buffer);
	model->update = _tree_handle_change;
	model->free = _tree_model_free;
	model->clone = _tree_model_clone;
	model->get_free_parameters = _tree_model_get_free_parameters;
	return model;
}

#include "treeio.h"
Model* new_TreeModel_from_json(json_node* node, Hashtable* hash){
	json_node* newick_node = get_json_node(node, "newick");
	json_node* file_node = get_json_node(node, "file");
	json_node* init_node = get_json_node(node, "init");
	json_node* patterns_node = get_json_node(node, "patterns");
	
	Model* mtree = NULL;
	
	if (newick_node != NULL) {
		char* newick = (char*)newick_node->value;
		Tree* tree = new_Tree(newick, true);
		char* id = get_json_node_value_string(node, "id");
		mtree = new_TreeModel(id, tree);
		Hashtable_add(hash, get_json_node_value_string(node, "parameters"), tree->distances);
		for (int i = 0; i < Parameters_count(tree->distances); i++) {
			Hashtable_add(hash, Parameters_name(tree->distances, i), Parameters_at(tree->distances, i));
		}
	}
	else if (file_node != NULL) {
        const char* filename = (const char*)file_node->value;
        char* tree_string = readTree(filename);
        Tree *tree = new_Tree( tree_string, true );
        free(tree_string);
        char* id = get_json_node_value_string(node, "id");
        mtree = new_TreeModel(id, tree);
        Hashtable_add(hash, get_json_node_value_string(node, "parameters"), tree->distances);
        for (int i = 0; i < Parameters_count(tree->distances); i++) {
            Hashtable_add(hash, Parameters_name(tree->distances, i), Parameters_at(tree->distances, i));
        }

	}
	else if (init_node != NULL) {
		json_node* algorithm_node = get_json_node(init_node, "algorithm");
		json_node* model_node = get_json_node(init_node, "model");
		json_node* patterns_node = get_json_node(init_node, "sitepattern");
		
		if (algorithm_node != NULL) {
			char* algorithm = (char*)algorithm_node->value;
		}
		char* patterns_ref = (char*)patterns_node->value;
		
		SitePattern* patterns = Hashtable_get(hash, patterns_ref+1);

		distancematrix_model tt = DISTANCE_MATRIX_UNCORRECTED;
		if (model_node != NULL) {
			const char* model = (char*)model_node->value;
			if (strcasecmp("raw", model) == 0) {
				tt = DISTANCE_MATRIX_UNCORRECTED;
			}
			else if (strcasecmp("jc69", model) == 0) {
				tt = DISTANCE_MATRIX_JC69;
			}
			else if (strcasecmp("k2p", model) == 0) {
				tt = DISTANCE_MATRIX_K2P;
			}
			else if (strcasecmp("kimura", model) == 0) {
				tt = DISTANCE_MATRIX_KIMURA;
			}
			else{
				exit(1);
			}
		}
		
		double** matrix = SitePattern_distance(patterns, tt);
		Tree* tree = new_NJ((const char**)patterns->names, patterns->size, matrix);
		json_node* id = get_json_node(node, "id");
		mtree = new_TreeModel((char*)id->value, tree);
		free_dmatrix(matrix, patterns->size);
		Hashtable_add(hash, get_json_node_value_string(node, "parameters"), tree->distances);
		for (int i = 0; i < Parameters_count(tree->distances); i++) {
			Hashtable_add(hash, Parameters_name(tree->distances, i), Parameters_at(tree->distances, i));
		}
	}
	else if (node->node_type != MJSON_STRING) {
		char* ref = (char*)node->value;
		mtree = Hashtable_get(hash, ref+1);
		mtree->ref_count++;
	}
	else{
		exit(1);
	}
	
	return mtree;
}

static void free_Tree_aux( Node *n ){
	if( n != NULL ){
		free_Tree_aux(n->left);
		free_Tree_aux(n->right);
		free_Node( n );
	}
}

void free_Tree( Tree *t){
	free_Tree_aux( t->root );
	if ( t->postorder != NULL ) free(t->postorder);
	if ( t->preorder != NULL ) free(t->preorder);
	t->postorder = NULL;
	t->preorder = NULL;
	free(t->nodes);
	free_Parameters(t->distances);
	free(t);
	t = NULL;
}

void clone_Tree_aux2( const Node *node, Node *newNode, Node *parent ){
	if ( node == NULL ) {
		return;
	}
	fprintf(stderr, "-- %s\n",node->name);
	newNode = clone_Node(node);
	/*newNode->parent = parent;
	 
	 if ( parent != NULL ){
	 if ( parent->left == NULL ) {
	 parent->left = newNode;
	 }
	 else if ( parent->right == NULL ){
	 parent->right = newNode;
	 }
	 else {
	 error("shousfdo\n");		
	 }
	 
	 fprintf(stderr, "%s %s\n",parent->name,newNode->name);
	 }*/
	
	clone_Tree_aux2(node->left, newNode->left, newNode);
	
	clone_Tree_aux2(node->right, newNode->right, newNode);
	
}

Node * clone_Tree_aux( Tree *tree, const Node *node, Node *parent ){
	if ( node == NULL ) {
		return NULL;
	}
	
	Node *newNode = clone_Node(node);
    tree->nodes[ newNode->id ] = newNode;
	newNode->parent = parent;
	newNode->left  = clone_Tree_aux(tree, node->left, newNode );
	newNode->right = clone_Tree_aux(tree, node->right, newNode );
	
	return newNode;
}

Tree * clone_Tree( const Tree *tree ){
	Tree *newTree = (Tree *)malloc(sizeof(Tree));
	assert(newTree);
	
	newTree->id = tree->id+1;
	
	newTree->root      = NULL;
	newTree->postorder = NULL;
	newTree->preorder  = NULL;
    newTree->nodes = (Node**)malloc( Tree_node_count(tree) * sizeof(Node*) );
    assert(newTree->nodes);
	
	newTree->nTips  = tree->nTips;
	newTree->nNodes = tree->nNodes;
	newTree->rooted = tree->rooted;
	
	
	newTree->root = clone_Tree_aux( newTree, tree->root, NULL);
	
	newTree->topology_changed = tree->topology_changed;
	
    Tree_update_topology(newTree);
	
	newTree->dated = tree->dated;
	newTree->time_mode = tree->time_mode;
	
	newTree->distances = new_Parameters(Tree_node_count(newTree)-2);
	Node* root = Tree_root(newTree);
	for (int i = 0; i < Tree_node_count(newTree); i++) {
		Node* node = Tree_node(newTree, i);
		if ( node != root && root->right != node) {
			Parameters_add(newTree->distances, node->distance);
		}
	}
	return newTree;
}

void _Tree_count_nodes( Node *node, int *tips, int *nodes){
	if ( node == NULL ) {
		return;
	}
	_Tree_count_nodes( Node_left(node), tips, nodes );
	_Tree_count_nodes( Node_right(node), tips, nodes );
	if ( Node_isleaf(node) ) {
		(*tips)++;
	}
	(*nodes)++;
}

Tree * clone_SubTree( const Tree *tree, Node *node ){
	Tree *newTree = (Tree *)malloc(sizeof(Tree));
	assert(newTree);
	
	newTree->id = tree->id+1;
	
	newTree->postorder = NULL;
	newTree->preorder  = NULL;
    newTree->nodes     = (Node**)malloc( Tree_node_count(tree) * sizeof(Node*) );
    assert(newTree->nodes);
	
	newTree->rooted = tree->rooted;
	
	newTree->root = clone_Tree_aux(newTree, node, NULL);
	newTree->nTips  = 0;
	newTree->nNodes = 0;
	_Tree_count_nodes( newTree->root, &(newTree->nTips), &(newTree->nNodes) );
	
    
	newTree->topology_changed = tree->topology_changed;
    Tree_update_topology(newTree);
    
	newTree->dated = tree->dated;
	
	newTree->distances = new_Parameters(Tree_node_count(newTree)-2);
	Node* root = Tree_root(newTree);
	for (int i = 0; i < Tree_node_count(newTree); i++) {
		Node* node = Tree_node(newTree, i);
		if ( node != root && root->right != node) {
			Parameters_add(newTree->distances, node->distance);
		}
	}
	
	return newTree;
}


char * Tree_stringify( Tree *tree ){
	StringBuffer *buffer = new_StringBuffer(1000);
	
	buffer = Tree_Bufferize( buffer, tree );
	
	char *final = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	//fprintf(stderr, "%d\n", (int)strlen(final));
	return final;
}

static StringBuffer * _to_nexus( StringBuffer *buffer, const Node *n );


StringBuffer * Tree_Bufferize( StringBuffer *buffer, Tree *tree ){
	StringBuffer_append_format(buffer, "(Tree:\n (id:\"tree%d\")\n(nexus:\"", tree->id);
	
	_to_nexus(buffer, tree->root );
	StringBuffer_append_string(buffer,"\")\n");
	
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	int i = 0;
	if ( nodes[0]->distance != NULL ) {
		StringBuffer_append_format(buffer,"(%s:\n",POSTFIX_DISTANCE);
		for ( i = 0; i < tree->nNodes; i++) {
			Parameter_bufferize(buffer, nodes[i]->distance);
			StringBuffer_append_char(buffer,'\n');
		}
		StringBuffer_append_string(buffer,")\n");
	}
	
	if ( nodes[0]->height != NULL ) {
		StringBuffer_append_format(buffer,"(%s:\n",POSTFIX_HEIGHT);
		for ( i = 0; i < tree->nNodes; i++) {
			Parameter_bufferize(buffer, nodes[i]->height);
			StringBuffer_append_char(buffer,'\n');
		}
		StringBuffer_append_string(buffer,")\n");
	}
	
	
	
	StringBuffer_append_char(buffer,')');
	
	return buffer;
}

void * Tree_SML_to_object( ObjectStore *store, SMLNode node ){
	fprintf(stderr, "Tree_SML_to_object ");
	Tree *tree = NULL;
	
	char *id = SML_get_data_of_child( node, "id");
	
	
	if ( (tree = ObjectStore_get_object(store, id )) != NULL ) {
		fprintf(stderr, " use Object store\n");
		return tree;
	}
	
	SMLNode tree_node = SML_get_element( node, "nexus");
	if ( tree_node != NULL ) {
		fprintf(stderr, " create object\n");
		tree = new_Tree( SML_get_data(tree_node), true );
		
		SMLNode distance_node = SML_get_element( node, POSTFIX_DISTANCE);
		if ( distance_node != NULL ) {
			Parameters *ps = Parameters_SML_to_object(distance_node);
			tree_init_distance_parameters( tree, ps );
		}
		
		SMLNode height_node = SML_get_element( node, POSTFIX_HEIGHT);
		if ( distance_node != NULL ) {
			Parameters *ps = Parameters_SML_to_object(height_node);
			tree_init_height_parameters( tree, ps );
		}
		
		char *pid = &id[0];
		while( *pid < '0' || *pid > '9' ){
			pid++;
		}
		tree->id = atoi( id );
		
		ObjectStore_add( store, id, tree);
	}
	
	
	return tree;
}

Node ** Tree_nodes( Tree *tree ){
    return tree->nodes;
}

Node * Tree_node( Tree *tree, int index ){
    return tree->nodes[index];
}

int Tree_node_count( const Tree *tree ){
	return tree->nNodes;
}

int Tree_tip_count( const Tree *tree ){
	return tree->nTips;
}


Node * Tree_root( Tree *tree ){
	return tree->root;
}

// use with caution
void Tree_set_root( Tree *tree, Node *root ){
    tree->root = root;
}

int Tree_id( const Tree *tree ){
	return tree->id;
}

void Tree_set_dated( Tree *tree, bool dated ){
	tree->dated = dated;
}

bool Tree_dated( const Tree *tree ){
	return tree->dated;
}

void Tree_set_id( Tree *tree, int id ){
	tree->id = id;
}

void Tree_post_to_preorder( Tree *tree, double *v ){
	double *a = dvector( Tree_node_count(tree) );
	Node **nodes = Tree_get_nodes(tree, PREORDER);
	for (int i = 0; i < Tree_node_count(tree); i++) {
		a[i] = v[ nodes[i]->postorder_idx ];
	}
	memcpy(v, a, Tree_node_count(tree)*sizeof(double) );
	free(a);
}


static void _rename_tree_postorder( Node *node, int *index ){
	if( node == NULL ) return;
	_rename_tree_postorder( node->left, index );
	_rename_tree_postorder( node->right, index );
	
	if( node->left != NULL ) Node_rename( node, *index );
	*index = *index + 1;
}

static void _rename_tree_preorder( Node *node, int *index ){
	if( node == NULL ) return;
	if( node->left != NULL ) Node_rename( node, *index );
	*index = *index + 1;
	_rename_tree_preorder( node->left, index );
	_rename_tree_preorder( node->right, index );
	
}

void rename_tree( Tree *tree, treeorder order ){
	int index = 0;
	switch ( order) {
		case POSTORDER:
			_rename_tree_postorder( tree->root, &index );
			break;
		case PREORDER:
			_rename_tree_preorder( tree->root, &index );
			break;
		default:
			break;
	}
}

StringBuffer * _to_nexus( StringBuffer *buffer, const Node *n ){
	if( n != NULL ){
		if( n->left != NULL ) buffer = StringBuffer_append_char(buffer, '(');
		else{
			buffer = StringBuffer_append_string(buffer, n->name);
			return buffer;
		}
		
		buffer = _to_nexus( buffer, n->left );
		buffer = StringBuffer_append_char(buffer, ',');
		_to_nexus( buffer, n->right );
		buffer = StringBuffer_append_char(buffer, ')');
	}
	return buffer;
}

static void _Tree_to_string_nexus_with_annotation_aux( StringBuffer *treebuffer, Tree *tree, const Node *n, int *count ){
    if( n == NULL ) return;
	if( !Node_isleaf(n) ) StringBuffer_append_char(treebuffer, '(');// fprintf(pf, "(");
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
            //fprintf(pf, "%d:[&%s]%f", ++(*count), buff->c, (n->parent->height->value - n->height->value) );
            StringBuffer_append_format(treebuffer,"%d:[&%s]%f", ++(*count), buff->c, (n->parent->height->value - n->height->value) );
            //fflush(pf);
            free_StringBuffer(buff);
        }
        else {
            //fprintf(pf, "%d:%f", ++(*count), (n->parent->height->value - n->height->value) );
            StringBuffer_append_format(treebuffer, "%d:%f", ++(*count), (n->parent->height->value - n->height->value) );
        }
		return;
	}
	
	_Tree_to_string_nexus_with_annotation_aux( treebuffer, tree, n->left, count );
	
	//fprintf(pf, ",");
    StringBuffer_append_char(treebuffer, ',');
	
	_Tree_to_string_nexus_with_annotation_aux( treebuffer, tree, n->right, count );
    
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
            //fprintf(pf, ")[&%s];", buff->c );
            StringBuffer_append_format(treebuffer, ")[&%s];", buff->c );
        }
        else {
            //fprintf(pf, "):[&%s]%f", buff->c, Node_time_elapsed((Node*)n) );
            StringBuffer_append_format(treebuffer, "):[&%s]%f", buff->c, Node_time_elapsed((Node*)n) );
        }
        free_StringBuffer(buff);
    }
    else {
        if( Node_isroot(n) ) StringBuffer_append_string(treebuffer, ");");// fprintf(pf, ");" );
        else StringBuffer_append_format(treebuffer, "):%f", Node_time_elapsed((Node*)n) );//fprintf(pf, "):%f", Node_time_elapsed((Node*)n) );
    }
}

char * Tree_to_string_nexus_with_annotation( Tree *tree ){
    int cunt = 0;
    StringBuffer *treebuffer = new_StringBuffer(100);
    _Tree_to_string_nexus_with_annotation_aux(treebuffer, tree, Tree_root(tree), &cunt);
    char *string = StringBuffer_tochar(treebuffer);
    free_StringBuffer(treebuffer);
    return string;
}

void Tree_StringBuffer_nexus_with_annotation( StringBuffer *treebuffer, Tree *tree ){
    int cunt = 0;
    _Tree_to_string_nexus_with_annotation_aux(treebuffer, tree, Tree_root(tree), &cunt);
}

// topology and distance
void Tree_newick_distance( StringBuffer *buffer, Node *n ){
	if( n != NULL ){
		if( !Node_isleaf(n) ) StringBuffer_append_char(buffer, '(');
		else{
            StringBuffer_append_format(buffer, "%s:%f", n->name, Node_distance(n));
		}
		
		Tree_newick_distance( buffer, Node_left(n) );
		StringBuffer_append_char(buffer, ',');
		Tree_newick_distance( buffer, Node_right(n) );
		StringBuffer_append_char(buffer, ')');
        if( !Node_isroot(n) ) StringBuffer_append_format(buffer, ":%f", n->distance->value);
	}
}

// topology only
void Tree_to_newick( StringBuffer *buffer, Node *n ){
	if( n != NULL ){
		if( !Node_isleaf(n) ) StringBuffer_append_char(buffer, '(');
		else{
			StringBuffer_append_string(buffer, n->name);
		}
		
		Tree_newick_distance( buffer, Node_left(n) );
		StringBuffer_append_char(buffer, ',');
		Tree_newick_distance( buffer, Node_right(n) );
		StringBuffer_append_char(buffer, ')');
	}
}

static void _Tree_init_depth_aux( Node *node ){
    if( node == NULL )return;
    
    if( !Node_isroot(node) ) node->depth =  node->parent->depth+1;
    
    _Tree_init_depth_aux( Node_left(node) );
    _Tree_init_depth_aux( Node_right(node) );
}

void Tree_init_depth( Tree *tree ){
	Tree_root(tree)->depth = 0;
	_Tree_init_depth_aux( Tree_root(tree) );
}

// What happens when memory is alraady allocated for a parameter
void tree_init_distance_parameters( Tree *t, Parameters *ps ){
	/*if ( t->parameters == NULL ) {
		t->parameters = new_Parameters( ps->count );
	}*/
	int i;
	Node ** nodes = Tree_get_nodes( t, POSTORDER );
	for ( i = 0; i < t->nNodes; i++ ) {
		nodes[i]->distance = Parameters_at(ps,i);
	}
}

void tree_init_height_parameters( Tree *t, Parameters *ps ){
	int i;
	Node ** nodes = Tree_get_nodes( t, POSTORDER );
	for ( i = 0; i < t->nNodes; i++ ) {
		nodes[i]->height = Parameters_at(ps,i);
	}
}

static void preorder( Node *n,  Node **nodes, bool tipsOnly, int *pos ){
	if( n == NULL ) return;
	if( (tipsOnly && Node_isleaf(n) ) || !tipsOnly ){
        //n->preorder_idx = *pos;
        nodes[(*pos)++] = n;
    }
	preorder( n->left,  nodes, tipsOnly, pos );
	preorder( n->right, nodes, tipsOnly, pos );
}

static void postorder( Node *n,  Node **nodes, bool tipsOnly, int *pos ){
	if( n == NULL ) return;
	postorder( n->left,  nodes, tipsOnly, pos );
	postorder( n->right, nodes, tipsOnly, pos );
	if( (tipsOnly && Node_isleaf(n) ) || !tipsOnly ){
        //n->postorder_idx = *pos;
        nodes[(*pos)++] = n;
    }
}

static void _postorder_generic_aux( Node *n, void *data, void (*action)( Node *node, void *data )){
	if( n == NULL ) return;
	_postorder_generic_aux( n->left, data, action);
	_postorder_generic_aux( n->right, data, action);
	(*action)(n,data);
}

void postorder_generic( Tree *tree, void *data, void (*action)( Node *node, void *data )){
	_postorder_generic_aux(tree->root, data, action);
}

static void _preorder_generic_aux( Node *n, void *data, void (*action)( Node *node, void *data )){
	if( n == NULL ) return;
	(*action)(n,data);
	_preorder_generic_aux( n->left, data, action);
	_preorder_generic_aux( n->right, data, action);
}

void pretorder_generic( Tree *tree, void *data, void (*action)( Node *node, void *data )){
	_preorder_generic_aux(tree->root, data, action);
}

void Tree_set_topology_changed( Tree *tree ){
    tree->topology_changed = true;
}

bool Tree_topology_changed( const Tree *tree ){
    return tree->topology_changed;
}

void _update_order_nodes( Tree *t ){
	int i;
	Node ** nodes = Tree_get_nodes( t, POSTORDER );
	for ( i = 0; i < t->nNodes; i++ ) {
		nodes[i]->postorder_idx = i;
	}
	
	nodes = Tree_get_nodes( t, PREORDER );
	for ( i = 0; i < t->nNodes; i++ ) {
		nodes[i]->preorder_idx = i;
	}
}

void Tree_update_topology( Tree *tree ){
    create_node_list( tree, POSTORDER );
    create_node_list( tree, PREORDER );
    tree->topology_changed = false;
    _update_order_nodes( tree );
    Tree_init_depth(tree);
}

void create_node_list( Tree *atree, treeorder order ){
	if ( order == POSTORDER && atree->postorder == NULL) {
		atree->postorder = (Node **)malloc( atree->nNodes * sizeof(Node*) );
		assert(atree->postorder);
	}
	else if ( order == PREORDER && atree->preorder == NULL) {
		atree->preorder = (Node **)malloc( atree->nNodes * sizeof(Node*) );
		assert(atree->preorder);
	}
	else if ( order == POSTORDER || order == PREORDER){
		
	}
	else {
		fprintf(stderr, "WTF create_node_list order:%d\n",order);
		exit(1);
	}
    
	
	int foo = 0;
	switch( order ){
		case PREORDER:
			preorder( atree->root,  atree->preorder, false, &foo );
			break;
		case POSTORDER:
			postorder( atree->root,  atree->postorder, false, &foo );
			break;
		default:
			error("create_node_list: not yet implemented\n");
	}
}

Node ** get_tips( Tree *atree, treeorder order ){
    if( atree->topology_changed ){
        Tree_update_topology(atree);
	}
    
	Node **tips = (Node **)malloc( atree->nTips * sizeof(Node*) );
	assert(tips);
	int foo = 0;
	switch( order ){
		case PREORDER:
			preorder( atree->root,  tips, true, &foo );
			return tips;
		case POSTORDER:
			postorder( atree->root,  tips, true, &foo );
			return tips;
		default:
			return NULL;
	}
}




Node ** Tree_get_nodes( Tree *atree, treeorder order ){
	if( atree->topology_changed ){
        Tree_update_topology(atree);
	}
	switch ( order ) {
		case PREORDER:
			return atree->preorder;
			break;
		case POSTORDER:
			return atree->postorder;
			break;
		default:
			break;
	}
	return NULL;
}

Node * Tree_get_node(Tree *atree, treeorder order, int index){
    if( atree->topology_changed ){
		Tree_update_topology(atree);
	}
    
	switch ( order ) {
		case PREORDER:
			return atree->preorder[index];
			break;
		case POSTORDER:
			return atree->postorder[index];
			break;
		default:
			break;
	}
    exit(2);
	return NULL;
}

Node * Tree_get_node_by_name(Tree *tree, const char *name ){
    
    Node **nodes = Tree_nodes(tree);
    for ( int i = 0; i< Tree_node_count(tree); i++ ) {
        if( strcmp(name, nodes[i]->name) == 0 ){
            return nodes[i];
        }
    }
    return NULL;
}

// swap their parents and update the tree structure
void Tree_rearrange( Tree *tree, Node *node1, Node *node2 ){
    Node_swap_parents(node1,node2);
    
	Tree_update_topology(tree);
}

// exchange nodes but do not update structure
void Tree_swap_parents_by_name( Tree *tree, char *name1, char *name2 ){
    Node *node1 = Tree_get_node_by_name(tree, name1);
    Node *node2 = Tree_get_node_by_name(tree, name2);
    
    Node_swap_parents(node1,node2);
}

int Tree_distance_nodes( Tree *tree, Node *node1, Node *node2 ){
    return Node_graph_distance(node1, node2);
}

#pragma mark -
// MARK: Distance methods


// assume no bounds
void Tree_scale_distance( Tree *tree, double scale ){
	Node **nodes = Tree_nodes( tree );
	for (int i = 0; i < Tree_node_count(tree); i++) {
        if( Node_isroot(nodes[i]) ) continue;
		Node_set_distance(nodes[i], Node_distance( nodes[i] )*scale);
	}
}

void Tree_branch_length_to_vector( Tree *tree, double *distances ){
	Node **nodes = Tree_nodes( tree );
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		distances[i] = Node_distance(nodes[i]);
	}
}


void Tree_vector_to_branch_length( Tree *tree, const double *distances ){
	Node **nodes = Tree_nodes( tree );
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		Node_set_distance(nodes[i], distances[i]);
	}
}

double Tree_distance_to_root( const Node *node ){
	Node *n = (Node*)node;
	double y = 0;
	while ( !Node_isroot(n) ) {
		y += Node_distance(n);
		n = Node_parent(n);
	}
	return y;
}

double Tree_distance_between_nodes( const Node *node, const Node *ancestor ){
	Node *n = (Node*)node;
	double y = 0;
	while ( n != ancestor ) {
		y += Node_distance(n);
		n = Node_parent(n);
	}
	return y;
}


double Tree_length( const Tree *tree ){
	double len = 0;
	Node **nodes = Tree_nodes((Tree *)tree);
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
        if( Node_isroot(nodes[i]) ) continue;
		len += Node_distance(nodes[i]);
	}
	return len;
}

// Scale tree length to total length
// Ignore constraints
double Tree_scale_total_length( Tree *tree, double total ){
	double len = 0;
	Node **nodes = Tree_nodes(tree);
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
        if( Node_isroot(nodes[i]) ) continue;
		len += Node_distance(nodes[i]);
	}
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
        if( Node_isroot(nodes[i]) ) continue;
		Node_set_distance(nodes[i], Node_distance(nodes[i])*total/len );
	}
	return len;
}


static void _Tree_copy_distances_aux( Node *src, Node *dst){
    if(src == NULL ) return;
    Node_set_distance( dst, Node_distance(src) );
    _Tree_copy_distances_aux(Node_left(src),Node_left(dst));
    _Tree_copy_distances_aux(Node_right(src),Node_right(dst));
}


// copy distances but does not check for bounds
void Tree_copy_distances( const Tree* src, Tree *dst ){
    _Tree_copy_distances_aux( Tree_root((Tree*)src), Tree_root(dst) );
	
}

#pragma mark -
// MARK: Height methods

void Tree_heights_to_vector( Tree *tree, double *heights ){
	Node **nodes = Tree_nodes( tree );
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		heights[i] = Node_height( nodes[i] );
	}
}

void Tree_vector_to_heights( const double *heights, Tree *tree ){
	Node **nodes = Tree_nodes( tree );
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		Node_set_height( nodes[i], heights[i] );
	}
}

// Initialize heights using branch lengths(subst) and rate
void Tree_init_heights_heterochronous( Tree *tree, const double rate, bool convert ){
	double max = 0;
	double youngest = -INFINITY; //makes sense only when convert==TRUE
    if ( !convert ) {
        youngest = INFINITY;
    }
	
	Node **nodes = Tree_get_nodes( tree, POSTORDER );
	
	double *distances = dvector(Tree_node_count(tree));
	
	Node *n = NULL;
	for (int i = 0; i < Tree_node_count(tree); i++) {
		n = nodes[i];
		double tot_bl = 0;
		while ( !Node_isroot(n) ) {
			tot_bl += Node_distance(n);
			n = Node_parent(n);
		}
		distances[i] = tot_bl;
		max = dmax(max, tot_bl);
		if ( Node_isleaf(nodes[i]) ){
            if ( convert ) {
                youngest = dmax(youngest, Node_time(nodes[i]));
            }
            else {
                youngest = dmin(youngest, Node_time(nodes[i]));
            }
			
		}
	}
	
	for (int i = 0; i < Tree_node_count(tree); i++) {
		if ( Node_isleaf(nodes[i]) ){
			if ( nodes[i]->height == NULL ) nodes[i]->height = new_Parameter_with_postfix(nodes[i]->name, POSTFIX_HEIGHT, nodes[i]->time, new_Constraint(nodes[i]->time, nodes[i]->time)); // constraint is temporary
			Node_set_height( nodes[i], nodes[i]->time );
			Parameter_set_estimate(nodes[i]->height, false);
			
			if ( convert ) {
				nodes[i]->height->value = youngest - nodes[i]->height->value;
				Parameter_set_bounds(nodes[i]->height, Parameter_value(nodes[i]->height), Parameter_value(nodes[i]->height) );
			}
            else {
				nodes[i]->height->value = nodes[i]->height->value - youngest;
				Parameter_set_bounds(nodes[i]->height, Parameter_value(nodes[i]->height), Parameter_value(nodes[i]->height) );                
            }
		}
		else {
			if ( nodes[i]->height == NULL ) nodes[i]->height = new_Parameter_with_postfix(nodes[i]->name, POSTFIX_HEIGHT, ((max - distances[i])/rate), new_Constraint(0, 0));
			Node_set_height(nodes[i],((max - distances[i])/rate));
			
			double max_children = dmax( Node_height(Node_left(nodes[i])), Node_height(Node_right(nodes[i])));
			if( nodes[i]->height->value < max_children ){
				nodes[i]->height->value = max_children + 0.001;
			}
		}

	}
	free(distances);
	
	Tree_constraint_heights(tree);
}

// initialize the bounds of the calibrated nodes
bool parse_calibrations( Tree *tree ){
	StringBuffer *buffer = new_StringBuffer(100);
	Node **nodes = Tree_get_nodes( tree, POSTORDER );
	bool isCal = false;
	
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {		
		
		if ( nodes[i]->info != NULL ) {
            //printf("%s %s\n", nodes[i]->name, nodes[i]->info);
			int pos = 0;
			pos = String_index_of_str(nodes[i]->info, CALIBRATION_TAG);
			if ( pos != -1 ) {
				isCal = true;
				//fprintf(stderr, "%s %s\n", nodes[i]->name, nodes[i]->info);
				char *ptr = nodes[i]->info+pos;
				while ( *ptr != '{' ) {ptr++;}
				ptr++;
				StringBuffer_empty(buffer);
				
				while ( *ptr != ',' ) {
					StringBuffer_append_char(buffer, *ptr);
					ptr++;
				}
				
				if ( buffer->length > 0 ) {
					double lower = atof(buffer->c);
					//fprintf(stderr,"lower %s %f\n", buffer->c, lower);
					Parameter_set_lower(nodes[i]->height, lower );
					Constraint_set_flower(nodes[i]->height->cnstr, lower);
					Constraint_set_lower_fixed(nodes[i]->height->cnstr, true);
				}
				
				StringBuffer_empty(buffer);
				
				ptr++;
				
				while ( *ptr != '}' ) {
					StringBuffer_append_char(buffer, *ptr);
					ptr++;
				}
				if ( buffer->length > 0 ) {
					double upper = atof(buffer->c);
					Parameter_set_upper(nodes[i]->height, upper);
					Constraint_set_fupper(nodes[i]->height->cnstr, upper);
					Constraint_set_upper_fixed(nodes[i]->height->cnstr, true);
				}
				
				
				
				//				if ( Parameter_upper(nodes[i]->height) == Parameter_lower(nodes[i]->height) ) {
				//					Parameter_set_fixed(nodes[i]->height, true);
				//				}
				if ( Constraint_lower_fixed(nodes[i]->height->cnstr) && Constraint_upper_fixed(nodes[i]->height->cnstr) && Parameter_upper(nodes[i]->height) == Parameter_lower(nodes[i]->height) ) {
					Parameter_set_value( nodes[i]->height, Parameter_upper(nodes[i]->height) );
					Parameter_set_estimate(nodes[i]->height, false);
					//fprintf(stderr, "%s %s\n", nodes[i]->name, nodes[i]->info);
				}
				else if ( Constraint_lower_fixed(nodes[i]->height->cnstr) && Constraint_upper_fixed(nodes[i]->height->cnstr) ) {
					Parameter_set_value(nodes[i]->height, (Parameter_upper(nodes[i]->height) - Parameter_lower(nodes[i]->height))*0.5 + Parameter_lower(nodes[i]->height));
				}
				else if( Constraint_lower_fixed(nodes[i]->height->cnstr) ){
					Parameter_set_value(nodes[i]->height, Parameter_lower(nodes[i]->height) );
				}
				else if( Constraint_upper_fixed(nodes[i]->height->cnstr) ){
					Parameter_set_value(nodes[i]->height, Parameter_upper(nodes[i]->height) );
					
				}
				
				
			}
		}
		
	}
	free_StringBuffer(buffer);
	return isCal;
}

// used Tree_init_heights_homochronous
static void _init_height_homochronous( const Node *thenode, Node *node, double max, double *distances, bool *set ){
	if ( node == NULL || Node_isleaf(node) ) {
		return;
	}
	
	Constraint *cnstr = node->height->cnstr;
	if ( thenode == node || !(Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr)) ) {
		_init_height_homochronous(thenode, node->left, max, distances, set);
		_init_height_homochronous(thenode, node->right, max, distances, set);
		
		if ( thenode != node ){
			//Parameter_set_value(node->height, (Node_height(thenode) * (max-distances[node->postorder_idx])/max) );
			double max_children = dmax( Node_height(Node_left(node)), Node_height(Node_right(node)));
            double max_children_d = dmax( Node_distance(Node_left(node)), Node_distance(Node_right(node)));
            
            Node_set_height(node, max_children+ max_children_d/max);
            
            //			if( node->height->value < max_children ){
            //				node->height->value = max_children + 0.001;
            //			}
            printf("-- %s %f\n", node->name, Node_height(node));
			set[node->postorder_idx] = true;
		}
		
	}
}

static void _init_height_homochronous2( const Node *thenode, Node *node, double rate, double maxd, bool *set ){
	if ( node == NULL || Node_isleaf(node) ) {
		return;
	}
	
	Constraint *cnstr = node->height->cnstr;
	if ( thenode == node || !(Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr)) ) {
		_init_height_homochronous2(thenode, node->left, rate, maxd, set);
		_init_height_homochronous2(thenode, node->right, rate, maxd, set);
		
		if ( thenode != node ){
			//Parameter_set_value(node->height, (Node_height(thenode) * (max-distances[node->postorder_idx])/max) );
			double max_children = dmax( Node_height(Node_left(node)), Node_height(Node_right(node)));
            double min_children_d = dmin( Node_distance(Node_left(node)), Node_distance(Node_right(node)));
            double scaler = Node_height(thenode) * rate/ maxd ;
            
            Node_set_height(node, max_children+ min_children_d*scaler/rate);
            
            //printf("-- %s %f scaler %f maxd %f max_children_d %f rate %f add %f\n", node->name, Node_height(node), scaler, maxd, min_children_d, rate,((min_children_d * scaler) /rate));
			set[node->postorder_idx] = true;
		}
		
	}
}



// Initialize heights using branch lengths(subst) and rate
void Tree_init_heights_homochronous( Tree *tree, const double rate ){

	Node **nodes = Tree_get_nodes( tree, POSTORDER );
	
	double *distances = dvector(Tree_node_count(tree)); // distance to root of each node
	
	bool *set = bvector(Tree_node_count(tree)); // if true then the height is already set because of a calibration
	
    double bl_tot  = Node_distance(Node_left(Tree_root(tree))) + Node_distance(Node_right(Tree_root(tree)));
    double bl_half = bl_tot/2;
    
    Node_set_distance(Node_left(Tree_root(tree)), bl_half);
    Node_set_distance(Node_right(Tree_root(tree)), bl_tot-bl_half);
    
	// find longest branch form root to tip
	double max_r = 0.0; // maximum distance to root

	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		if ( nodes[i]->height == NULL ) {
			nodes[i]->height = new_Parameter_with_postfix(nodes[i]->name, POSTFIX_HEIGHT, 0.0, new_Constraint(0.0, 0.0));
		}
        else if(Node_isleaf(nodes[i]) ) {
            Node_set_height(nodes[i], 0);
			Parameter_set_estimate(nodes[i]->height, false);
        }
		
		Node *n = nodes[i];
		double tot_bl = 0.0;
		while ( !Node_isroot( n ) ) {
			tot_bl += Node_distance(n);
			n = Node_parent(n);
		}
		max_r = dmax( tot_bl, max_r );
		distances[i] = tot_bl;
	}
	
	
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {		
		if ( Constraint_lower_fixed(nodes[i]->height->cnstr) && Constraint_upper_fixed(nodes[i]->height->cnstr) ) {
			double max = 0;
			for ( int j = i-1; j >= 0; j-- ) {
				if ( Node_isancestor( nodes[j], nodes[i]) ) {
					max = dmax( distances[j], max );
				}
				else break;
			}

			set[i] = true;
			
			//_init_height_homochronous(nodes[i], nodes[i], rate, distances, set);
            _init_height_homochronous2(nodes[i], nodes[i], rate,max, set);
			
			continue;
			
			for ( int j = i-1; j >= 0; j-- ) {
				if ( Node_isleaf(nodes[j]) || set[j] ) continue;
				if ( Node_isancestor( nodes[j], nodes[i]) ) {
					Parameter_set_value(nodes[j]->height, (Node_height(nodes[i]) * (max-distances[j])/max) );
					double max_children = dmax( Node_height(Node_left(nodes[j])), Node_height(Node_right(nodes[j])));
					if( nodes[j]->height->value < max_children ){
						nodes[j]->height->value = max_children + 0.001;
					}
					set[j] = true;
				}
				else break;
			}
			
			
			
		}
		
	}
	free(distances);

	
	// Set up all the other nodes up to the root
	if ( !set[ Tree_node_count(tree)-1 ] ) {
//		double max_height = 0.0;
//		double index_max = 0;
//		for (int i = 0; i < Tree_node_count(tree)-1; i++) {
//			if ( Node_height(nodes[i]) > max_height ) {
//				index_max = i;
//				max_height = Node_height(nodes[i]);
//			}
//		}
		
		for (int i = 0; i < Tree_node_count(tree); i++) {
			if ( set[i] || Node_isleaf(nodes[i]) ) continue;
//			Parameter_set_value(nodes[i]->height, ((max_r - d_r[i])/rate) );
//			//Parameter_set_value(nodes[j]->height, (Node_height(nodes[i]) * distances[i]/max_r) );
//			double max_children = dmax( Node_height(Node_left(nodes[i])), Node_height(Node_right(nodes[i])));
//			if( nodes[i]->height->value < max_children ){
//				nodes[i]->height->value = max_children + 0.001;
//			}
			double h = 0.0;
			if ( Node_height(Node_left(nodes[i])) > Node_height(Node_right(nodes[i]))) {
				h = Node_height(Node_left(nodes[i])) + Node_distance(Node_left(nodes[i]))/rate;
			}
			else {
				h = Node_height(Node_right(nodes[i])) + Node_distance(Node_right(nodes[i]))/rate;
			}
			Parameter_set_value(nodes[i]->height, h );
			
		}
	}
	free(set);
	
	Tree_constraint_heights(tree);
}

void summarize_calibrations( Tree *tree ){
	double rate = 0.0;
	int count = 0;
	double max = 0.0;
	double min = HUGE_VALF;
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		if ( !Node_isleaf( nodes[i] ) ) {
			Constraint *cnstr = nodes[i]->height->cnstr;
			if ( Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr) ) {
				
				fprintf(stderr, "\n%s %s\n", Node_name(nodes[i]), nodes[i]->info);
				for ( int j = i-1; j >= 0; j-- ) {
					if ( !Node_isleaf(nodes[j]) ) continue;
					if ( !Node_isancestor( nodes[j], nodes[i]) ) break;
					Node *node = nodes[j];
					fprintf(stderr, "\n%s \n", Node_name(node) );
					double d = 0.0;
					while ( node != nodes[i] ) {
						d += Node_distance(node);
						node = Node_parent(node);
					}
					fprintf(stderr, "d = %f\n", d);
					if ( Constraint_lower_fixed(cnstr) ) {
						double r = d/Constraint_flower(cnstr);
						fprintf(stderr, "Lower rate = %f\n", r );
						rate += r;
						count++;
						max = dmax(max, r);
					}
					if ( Constraint_upper_fixed(cnstr) ) {
						double r = (d/Constraint_fupper(cnstr));
						fprintf(stderr, "Upper rate = %f\n", r );
						rate += r;
						count++;
						min = dmin(min, r);
					}
					
				}
				
			}
		}
	}
	fprintf(stderr, "\nMin rate:        %f\n",min);
	fprintf(stderr, "Mean rate: %f\n",(rate/count));
	fprintf(stderr, "Max rate:        %f\n",max);
}

bool Tree_is_calibrated( Tree *tree ){
	Node **nodes = Tree_nodes(tree);
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		if ( !Node_isleaf( nodes[i] ) ) {
			Constraint *cnstr = nodes[i]->height->cnstr;
			if ( Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr) ) {
				return true;
			}
		}
	}
	return false;
}

double mean_rate_calibrations( Tree *tree ){
	double rate = 0.0;
	int count = 0;
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		if ( !Node_isleaf( nodes[i] ) ) {
			Constraint *cnstr = nodes[i]->height->cnstr;
			if ( Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr) ) {
				
				for ( int j = i-1; j >= 0; j-- ) {
					if ( !Node_isleaf(nodes[j]) ) continue;
					if ( !Node_isancestor( nodes[j], nodes[i]) ) break;
					Node *node = nodes[j];
					
					double d = 0.0;
					while ( node != nodes[i] ) {
						d += Node_distance(node);
						node = Node_parent(node);
					}
					
					if ( Constraint_lower_fixed(cnstr) ) {
						double r = d/Constraint_flower(cnstr);
						rate += r;
						count++;
					}
					if ( Constraint_upper_fixed(cnstr) ) {
						double r = (d/Constraint_fupper(cnstr));
						rate += r;
						count++;
					}
					
				}
				
			}
		}
	}
	return (rate/count);
}


// zero is the present, 10 is the past
bool parse_dates( Tree *tree ){
	char *pch = NULL;
	bool isDated = true;
	Node **nodes = Tree_get_nodes( tree, POSTORDER );
	
	for ( int i = 0; i < tree->nNodes; i++ ) {
		if( Node_isleaf(nodes[i]) ){
            int len = strlen(nodes[i]->name)-1;
			int j = len;
			pch = &nodes[i]->name[j];
			while ( *pch == '.' || (*pch >= 48 && *pch <= 57) ) {
				//if ( j == 0 || *pch == '_' || *pch == '/' ) {
                if ( j == 0 ) {
					break;
				}
				pch--;
				j--;
			}
			//if ( j == strlen(nodes[i]->name)-1 || ( j != strlen(nodes[i]->name)-1 && *pch != '_') ) {
            if ( j == len ) {
				isDated = false;
				break;
			}		
		}		
	}
	
	if ( isDated ) {
		StringBuffer *buffer = new_StringBuffer(10);
		for ( int i = 0; i < tree->nNodes; i++ ) {
			if( Node_isleaf(nodes[i]) ){
				int j = strlen(nodes[i]->name)-1;
				pch = &nodes[i]->name[j];
				StringBuffer_empty(buffer);
				//while ( *pch != '_' ) pch--;
                
                while ( *pch == '.' || (*pch >= 48 && *pch <= 57) ) pch--;
                pch++;
                
				while ( *pch != '\0' ) {
					StringBuffer_append_char(buffer, *pch);
					pch++;					
				}
				nodes[i]->time = atof( buffer->c);
				//fprintf(stderr, "%s %f %s\n", nodes[i]->name, nodes[i]->time, buffer->c);
			}
		}
		free_StringBuffer(buffer);
	}

	return isDated;
}



void Tree_constrain_height( Node *node ){
	Constraint *cnstr = node->height->cnstr;
	if ( Parameter_estimate(node->height) == false ) {
		return;
	}
	
	// Always assume that tips have a date
	if( Node_isleaf( node ) ){
		Constraint_set_bounds(cnstr, Parameter_value(node->height), Parameter_value(node->height) );
		Parameter_set_estimate(node->height, false);
	}
	else {
		double lower = dmax(Parameter_value( Node_left(node)->height ), Parameter_value( Node_right(node)->height ) );
		if( Constraint_lower_fixed(cnstr) ){
			lower = dmax( lower, Constraint_flower(cnstr) );
		}
		
		//double upper = ( Node_isroot(node) ? HEIGHT_MAX : Parameter_value( Node_parent(node)->height) );
		double upper = ( Node_isroot(node) ? Parameter_value(node->height)*2.0 : Parameter_value( Node_parent(node)->height) );
		if( Constraint_upper_fixed(cnstr) ){
			upper = dmin( upper, Constraint_fupper(cnstr) );
		}
		
		
		Constraint_set_lower(cnstr, lower);
		Constraint_set_upper(cnstr, upper );		
	}
}

void Tree_constrain_neighbor_heights( Node *node ){
	if( node->parent != NULL ) Tree_constrain_height(node->parent);
	if( node->left != NULL ) Tree_constrain_height(node->left);
	if( node->right != NULL ) Tree_constrain_height(node->right);
}


static void _Tree_constraint_heights_aux( Node *node ){
    if( node == NULL ) return;
    _Tree_constraint_heights_aux( Node_left(node));
    _Tree_constraint_heights_aux( Node_right(node));
    Tree_constrain_height( node );
}

// Done in postorder
void Tree_constraint_heights( Tree *tree ){
    _Tree_constraint_heights_aux(Tree_root(tree));
}



static void _Tree_check_constraint_heights_aux( Node *node ){
    if( node == NULL ) return;
    
        
    if( !Node_isleaf(node) ){
        if( ( Node_height(Node_left(node)) > Constraint_lower(node->height->cnstr) || Node_height(Node_right(node)) > Constraint_lower(node->height->cnstr) ) ){
            fprintf(stderr, "_Tree_check_constraint_heights_aux Lower Contraint %s\n", node->name);
            exit(1);
        }
        
        if( Node_height(Node_left(node)) > Node_height(node) || Node_height(Node_left(node)) > Node_height(node) ){
            fprintf(stderr, "_Tree_check_constraint_heights_aux Parameter %s\n", node->name);
            exit(1);
        }
    }
    
    if( !Node_isroot(node) ){
        if ( Node_height(Node_parent(node)) < Constraint_upper(node->height->cnstr) ){
            fprintf(stderr, "_Tree_check_constraint_heights_aux Upper Contraint %s\n", node->name);
            exit(1);
        }
        if( Node_height(Node_parent(node)) < Node_height(node)){
            fprintf(stderr, " _Tree_check_constraint_heights_aux Parameter parent %s\n", node->name);
            exit(1);
        }
    }
    _Tree_check_constraint_heights_aux( Node_left(node));
    _Tree_check_constraint_heights_aux( Node_right(node));
    
}


void Tree_check_constraint_heights( Tree *tree ){
    _Tree_check_constraint_heights_aux(Tree_root(tree));
}

void update_constraint_height_neighbors( Node *node ){
	if( node->parent != NULL){
		Parameter_set_lower(node->parent->height, dmax(node->parent->left->height->value, Parameter_value( Node_right(Node_parent(node))->height) ) );
	}
	if( !Node_isleaf(Node_left(node)) ) Parameter_set_upper( Node_left(node)->height, Parameter_value(node->height) );
	if( !Node_isleaf(Node_right(node)) ) Parameter_set_upper( Node_right(node)->height, Parameter_value(node->height) );
}

// times contains the distance between 2 nodes (e.g. years)
// The vector time is indexed by node ids
void Tree_heights_to_times( Node *node, double *times ){
	if( node == NULL ) return;
	Tree_heights_to_times( Node_left(node), times );
	Tree_heights_to_times( Node_right(node), times );
	
	if( !Node_isroot(node) ){
		times[ Node_id(node) ] = Node_height( Node_parent(node) ) - Node_height(node);
	}
	else {
		times[ Node_id(node) ] = 0.0;
	}
}


void Tree_print_parameters( Tree *tree ){
	Node **nodes = Tree_nodes( tree );
	for ( int i = 0; i < tree->nNodes; i++ ) {
		fprintf(stderr, "%d] postidx %d - %s height=%f [%f-%f] %d\n",
				i,nodes[i]->id, nodes[i]->name, nodes[i]->distance->value,
				Parameter_lower(nodes[i]->distance) ,
				Parameter_upper(nodes[i]->distance),
				Parameter_estimate(nodes[i]->distance) );
	}
}

// copy heights and constraints
void Tree_copy_height_constraints( Tree *src, Tree *dist ){
	Node **dist_nodes = Tree_nodes(dist);
	Node **src_nodes = Tree_nodes(src);
	for ( int i = 0; i < Tree_node_count(dist); i++ ) {
		//Constraint_set_lower( dist_nodes[i]->height->cnstr, Constraint_lower( src_nodes[i]->height->cnstr) );
		//Constraint_set_upper( dist_nodes[i]->height->cnstr, Constraint_upper( src_nodes[i]->height->cnstr) );
		
		Constraint_set_flower( dist_nodes[i]->height->cnstr, Constraint_flower( src_nodes[i]->height->cnstr) );
		Constraint_set_fupper( dist_nodes[i]->height->cnstr, Constraint_fupper( src_nodes[i]->height->cnstr) );
		
		Constraint_set_lower_fixed(dist_nodes[i]->height->cnstr, Constraint_lower_fixed(src_nodes[i]->height->cnstr));
		Constraint_set_upper_fixed(dist_nodes[i]->height->cnstr, Constraint_upper_fixed(src_nodes[i]->height->cnstr));
		Parameter_set_estimate(dist_nodes[i]->height, Parameter_estimate(src_nodes[i]->height));
	}
}

void Tree_copy_heights( Tree *src, Tree *dist ){
	Node **dist_nodes = Tree_nodes(dist);
	Node **src_nodes = Tree_nodes(src);
	for ( int i = 0; i < Tree_node_count(dist); i++ ) {
        Node_set_height(dist_nodes[i], Node_height(src_nodes[i]));
	}
}


void Tree_scale_heights( Node *node, const double scaler ){
	if ( node == NULL || Node_isleaf(node) ) {
		return;
	}
	
	Constraint *cnstr = node->height->cnstr;
	if ( Node_isroot(node) || !(Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr)) ) {
		Tree_scale_heights(node->left, scaler);
		Tree_scale_heights(node->right, scaler);
		
		if ( scaler > 1) {
			Node_set_height(node, Node_height( node )*scaler);
		}
		else if( scaler < 1 ){
			double height = Node_height( node );
			double max_height_son = dmax( Node_height( Node_left(node)), Node_height( Node_right(node) ) );
			double bl = (height  - max_height_son)*scaler + max_height_son;
			Node_set_height(node, bl);
		}
	}
}



#pragma mark -
#pragma mark Misc

//FIXME: should also update height parameters...
void Tree_reroot( Tree *tree, Node *node ){
	// node is already the root
	if ( Node_isroot(node)) {
		return;
	}
	// node is just below the root, just choose the middle of the branch
	if( Node_isroot( Node_parent(node) ) ){
		double midpoint = Node_distance(node)/2;
		Node_set_distance(node, Node_distance(node)-midpoint);
		
		Node *sibling = ( node == node->parent->left ? node->parent->right : node->parent->left);
		Node_set_distance(sibling, Node_distance(sibling)+midpoint);
	}
	else {
		
		Node *newroot = new_Node(NULL, "node0", 0);
		
		Node_addChild(newroot, node);
		Node_addChild(newroot, node->parent);
		
		double branchLength = Node_distance(node->parent); // save branch length
		double midpoint = Node_distance(node)/2;
		
		Node_set_distance(node, midpoint);
		Node_set_distance(node->parent, midpoint);
		
		Node *n = node->parent;     // the node to which we need to add its parent as a child recurssively
		Node *nparent = n->parent;  // the parent which is a ref to the rest of the old tree
		
		Node_removeChild(nparent, n);
		Node_removeChild(n, node);
		
		node->parent = newroot;
		n->parent = newroot;
        
		while ( !Node_isroot(nparent) ) {
			Node_addChild(n, nparent);
			
			Node *temp = nparent->parent;
			double bl = Node_distance(nparent);
			Node_set_distance(nparent, branchLength);
			branchLength = bl;
			
			nparent->parent = n;
			
			n = nparent;
			nparent = temp;
			
			Node_removeChild(nparent, n);
		}
		
		if( nparent->left == NULL){
			nparent = nparent->right;
		}
		else {
			nparent = nparent->left;
		}
		Node_set_distance(nparent, Node_distance(nparent) + branchLength );
		Node_addChild(n, nparent);
		nparent->parent = n;
		
        /*newroot->id = tree->root->id;
        
		free_Node(tree->root);
		tree->root = newroot;*/

        // new
        // The old root is reused
        Node *root = Tree_root(tree);
        Node_set_distance(root, -1);
        root->left = newroot->left;
        root->right = newroot->right;
        Node_set_parent(root->left, root);
        Node_set_parent(root->right, root);
        free_Node(newroot);

	}
	
    Tree_set_topology_changed(tree);
}





void compare_tree_aux( const Node *n1, const Node *n2){
	
	if( n1 == NULL && n2 == NULL){
		return;
	}
	else if ( n1 != NULL && n2 != NULL){}
	else {
		error("n1 and n2 !=");
	}
	assert(n1 != n2 );
	assert( strcmp(n1->name, n2->name) == 0 );
	
	if( n1->parent != NULL && n2->parent != NULL ) assert( strcmp(n1->name, n2->name) == 0 );
	else assert( n1->parent == NULL && n2->parent == NULL );
	
	assert( n1->id == n2->id );
	assert( n1->depth == n2->depth);
	assert( n1->time == n2->time);
	assert( n1->postorder_idx == n2->postorder_idx);
	assert( n1->preorder_idx == n2->preorder_idx);
	
	compare_parameter( n1->height, n2->height );
	compare_parameter( n1->distance, n2->distance );
		
	compare_tree_aux( n1->left, n2->left );
	compare_tree_aux( n1->right, n2->right );
}

void compare_tree( const Tree *t1, const Tree *t2 ){
	assert( t1->id == t2->id );
	assert( t1->nTips == t2->nTips );
	assert( t1->nNodes == t2->nNodes );
	assert( t1->rooted == t2->rooted );
	assert( t1->topology_changed == t2->topology_changed );
	assert( strcmp(t1->root->name, t2->root->name) == 0 );
	compare_tree_aux( t1->root, t2->root );
}

void print_compare_tree(  Tree *t1, Tree *t2){
	printf("ID %d %d\nnTips %d %d\nnNodes %d %d\nrooted %d %d\ntopology_changed %d %d\n",t1->id, t2->id, t1->nTips, t2->nTips, t1->nNodes, t2->nNodes, t1->rooted, t2->rooted,t1->topology_changed, t2->topology_changed);
	printf("root %s %s\n",t1->root->name, t2->root->name);
	Node **n1 = Tree_get_nodes( t1, POSTORDER);
	Node **n2 = Tree_get_nodes( t2, POSTORDER);
	for (int i = 0; i < t1->nNodes; i++) {
		printf("%s %s %d %d depth %d %d time %f %f post %d %d pre %d %d\n",
			   n1[i]->name,n2[i]->name,
			   n1[i]->id,n2[i]->id,
			   n1[i]->depth, n2[i]->depth,
			   n1[i]->time, n2[i]->time,
			   n1[i]->postorder_idx, n2[i]->postorder_idx,
			   n1[i]->preorder_idx, n2[i]->preorder_idx);
		
		Parameter_print( n1[i]->distance );
		Parameter_print( n2[i]->distance );
		
		Parameter_print( n1[i]->height );
		Parameter_print( n2[i]->height );
	}
}







