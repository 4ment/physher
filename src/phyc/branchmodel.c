/*
 *  branchmodel.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 12/8/10.
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

#include "branchmodel.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "tree.h"
#include "matrix.h"
#include "utils.h"
#include "mstring.h"
#include "parser.h"
#include "random.h"
#include "statistics.h"

#include "treeio.h"
#include "mathconstant.h"
#include "exponential.h"
#include "lognormal.h"
#include "parameters.h"
#include "node.h"




#define STRICT_POSTFIX   "strictclock"
#define LOCAL_POSTFIX    "localclock"
#define DISCRETE_POSTFIX "discreteclock"
#define RELAXED_POSTFIX  "relaxedclock"


static void _free_BranchModel( BranchModel *bm, bool remove_tree  ){
	if( bm->indicators != NULL ) free(bm->indicators);
	if( bm->map != NULL ) bm->map->free(bm->map);
	if( bm->unscaled_rates != NULL ) free(bm->unscaled_rates);
	if( bm->rates != NULL ) free_Parameters(bm->rates);
	if ( remove_tree && bm->tree != NULL ) {
		free_Tree(bm->tree);
	}
	free(bm);
}

/**************************** BRANCH MODEL ****************************/
// MARK: Public functions

BranchModel * BranchModel_init( Tree *tree, branchmodel name ){
	
	BranchModel *bm = (BranchModel *)malloc(sizeof(BranchModel));
	assert(bm);
	bm->id = 0;
	bm->tree = tree;
	bm->name = name;

	bm->rates = NULL;
	
	bm->get = NULL;
	bm->set = NULL;
    
	bm->free = _free_BranchModel;
	
	bm->need_update = true;
	
	// LOCAL clock
	bm->indicators = NULL;
	bm->unscaled_rates= NULL;
	bm->scalefactor = 0.0;
	
	// RELAXED
	bm->type = -1;
	
	// local and discrete
	bm->map = NULL;
	
	return bm;
}

BranchModel * new_BranchModel( Tree *tree, branchmodel type ){
	
	switch (type) {
		case CLOCK_STRICT:
			return new_StrictClock( tree );
			break;
		case CLOCK_LOCAL:
			return new_LocalClock( tree, 1 );
			break;
		case CLOCK_DISCRETE:
			return new_DiscreteClock( tree, 2 );
			break;
		default:
            break;
	}
	return new_StrictClock( tree );
}


#ifdef LISTENERS
static void _branchmodel_handle_change( Model *self, Model *model, int index ){
	// from the index we can try to tell treelikelihood to update the right part of the tree
	BranchModel *bm = (BranchModel*)self->obj;
	if(bm->name == CLOCK_STRICT){
		self->listeners->fire(self->listeners, self, -1 );
	}
	else if(bm->name == CLOCK_LOCAL){
		self->listeners->fire(self->listeners, self, index );
	}
	else{
		self->listeners->fire(self->listeners, self, -1 );
	}
}

static void _branch_model_free( Model *self ){
	if(self->ref_count == 1){
		printf("Free branch model %s\n", self->name);
		BranchModel *bm = (BranchModel*)self->obj;
		Model* mm = (Model*)self->data;
		mm->free(mm); // tree model
		
		//bm->free(bm, false);
		if( bm->indicators != NULL ) free(bm->indicators);
		bm->map->free(bm->map);
		if( bm->unscaled_rates != NULL ) free(bm->unscaled_rates);
		free_Parameters(bm->rates);
		free(bm);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _branch_model_clone( Model* self, Hashtable *hash ){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	Model* mtree = (Model*)self->data;
	Model* mtreeclone = NULL;
	// Tree may have been parsed already
	if (Hashtable_exists(hash, mtree->name)) {
		mtreeclone = Hashtable_get(hash, mtree->name);
	}
	else{
		mtreeclone = mtree->clone(mtree, hash);
		Hashtable_add(hash, mtree->name, mtreeclone);
	}
	BranchModel*bmclone = clone_BranchModel((BranchModel*)self->obj, (Tree*)mtreeclone->obj);
	for (int i = 0; i < Parameters_count(bmclone->rates); i++) {
		Hashtable_add(hash, Parameters_name(bmclone->rates, i), Parameters_at(bmclone->rates, i));
	}
	if(bmclone->map != NULL){
		Hashtable_add(hash, bmclone->map->name, bmclone->map);
	}
	Model* clone = new_BranchModel2(self->name, bmclone, mtreeclone);
	Hashtable_add(hash, clone->name, clone);
	return clone;
}

// BranchModel2 listen to the rate parameters
Model * new_BranchModel2( const char* name, BranchModel *bm, Model* tree){
	Model *model = new_Model(name, bm);
	for ( int i = 0; i < Parameters_count(bm->rates); i++ ) {
		Parameters_at(bm->rates, i)->listeners->add( Parameters_at(bm->rates, i)->listeners, model );
	}
	if(bm->map != NULL){
		bm->map->listeners->add(bm->map->listeners, model);
	}
	model->update = _branchmodel_handle_change;
	model->free = _branch_model_free;
	model->clone = _branch_model_clone;
	model->data = tree;
	tree->ref_count++;
	return model;
}
#endif


BranchModel * clone_BranchModel(const BranchModel *bm, Tree *tree ){
	BranchModel *newbm = (BranchModel *)malloc(sizeof(BranchModel));
	assert(newbm);
	
	newbm->name = bm->name;
	newbm->id   = bm->id;
	
	newbm->get  = bm->get;
	newbm->set  = bm->set;
    
	newbm->free = bm->free;
	
	newbm->indicators = NULL;
	newbm->map = NULL;
	newbm->unscaled_rates = NULL;
	newbm->rates = NULL;
	
	newbm->tree = tree;
	
	
	newbm->need_update = bm->need_update;
	
	newbm->rates = clone_Parameters(bm->rates);
	
	newbm->scalefactor = bm->scalefactor;
	
	if( bm->indicators != NULL ){
		newbm->indicators     = clone_bvector( bm->indicators, Tree_node_count(tree) );
	}
	
	if (bm->map != NULL ) {
		newbm->map = bm->map->clone(bm->map);
	}
    if( bm->unscaled_rates != NULL ){
		newbm->unscaled_rates = clone_dvector( bm->unscaled_rates, Tree_node_count(tree) );
    }
	newbm->type = bm->type;
	
	return newbm;
}

char *BranchModel_stringify( BranchModel *bm ){
	StringBuffer *buffer = new_StringBuffer(500);

	buffer = BranchModel_bufferize( buffer, bm );	
	
	char *final = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return final;
}

StringBuffer * BranchModel_bufferize( StringBuffer *buffer, BranchModel *bm ){
	
	buffer = StringBuffer_append_format(buffer, "(BranchModel:\n(id:\"bm%d\")\n(type: \"%d\")\n(Tree: &tree%d)\n",bm->id,bm->name, Tree_id(bm->tree));
	
	buffer = Parameters_SML_bufferize(buffer, bm->rates);
	
	buffer = StringBuffer_append_char(buffer,')');
	//if( bm->name == LOCAL) buffer = StringBuffer_append_format(buffer, "(local: %d)\n", bm->rates->count);
	
	return buffer;
}

void * BranchModel_SML_to_object( ObjectStore *store, SMLNode node ){
	fprintf(stderr, "BranchModel_SML_to_object\n");
	BranchModel *bm = NULL;
    Parameters *rates = NULL;
	
	char *tree_ref = SML_get_data_of_child( node, "Tree");
	
	Tree * tree = ObjectStore_get( store, tree_ref );
	
	char *type     = SML_get_data_of_child( node, "type");
	int itype = atoi(type);
	
	SMLNode node_p = NULL;
    switch (itype) {
        case CLOCK_STRICT:
            node_p = SML_get_element( node, STRICT_POSTFIX);
            rates = Parameters_SML_to_object( node_p );
            bm = new_StrictClock_with_parameters( tree,rates );
            break;
        case CLOCK_LOCAL:
            node_p = SML_get_element( node, LOCAL_POSTFIX);
            rates = Parameters_SML_to_object( node_p );
            bm = new_LocalClock_with_parameters( tree, rates );
            break;
        default:
            assert(0);
    }
	
	char *id = SML_get_data_of_child( node, "id");
	if ( id != NULL ) {
		while( *id < '0' || *id > '9' ){
			id++;
		}
		bm->id = atoi( id );
	}
	else bm->id = 0;
	
	return bm;
}

// Parameters are cloned, including the constraint
Parameters * BranchModel_save_rates( BranchModel *bm ){
	Parameters *rates = new_Parameters( Parameters_count(bm->rates) );
	for ( int i = 0; i < Parameters_count(bm->rates); i++ ) {
		Parameters_add(rates, clone_Parameter( Parameters_at(bm->rates, i)) );
	}
	return rates;
}

void BranchModel_restore_rates( BranchModel *bm, const Parameters *rates ){
	for ( int i = 0; i < Parameters_count(bm->rates); i++ ) {
		Parameters_set_value(bm->rates, i, Parameters_value(rates, i) );
		Parameters_set_lower(bm->rates, i, Parameters_lower(rates, i) );
		Parameters_set_upper(bm->rates, i, Parameters_upper(rates, i) );
		Parameters_set_estimate(bm->rates, i, Parameters_estimate(rates, i) );
	}
}

void BranchModel_rates_to_vector( BranchModel *bm, double *rates ){
	for ( int i = 0; i < Parameters_count(bm->rates); i++ ) {
		rates[i] = Parameters_value(bm->rates, i);
	}

}

void BranchModel_vector_to_rates( BranchModel *bm, const double *rates ){
	for ( int i = 0; i < Parameters_count(bm->rates); i++ ) {
		Parameters_set_value(bm->rates, i, rates[i]);
	}
}

void BranchModel_value_to_rates( const double rate, BranchModel *bm ){
	for ( int i = 0; i < Parameters_count(bm->rates); i++ ) {
		Parameters_set_value(bm->rates, i, rate);
	}
}

int BranchModel_n_rate( BranchModel *bm ){
    return Parameters_count(bm->rates);
}

#pragma mark -
// MARK: Strict clock

/**************************** STRICT CLOCK ****************************/

static void _set_strickclock( BranchModel *bm, const int index, const double value );
static double _get_strict_clock( BranchModel *bm , Node *node );

BranchModel * new_StrictClock( Tree *tree ){
	Parameters *rates = new_Parameters(1);
	Parameters_move(rates, new_Parameter_with_postfix("rate", STRICT_POSTFIX, 1., new_Constraint(BRANCHMODEL_RATE_MIN, BRANCHMODEL_RATE_MAX)));
	return new_StrictClock_with_parameters(tree,rates);
}

BranchModel * new_StrictClock_with_parameters( Tree *tree, const Parameters *rates ){
	BranchModel *bm = BranchModel_init(tree, CLOCK_STRICT);
	
	bm->get = _get_strict_clock;
	bm->set = _set_strickclock;
	
	bm->rates = new_Parameters(1);
	Parameters_add_parameters(bm->rates, rates);
	
	bm->need_update = false; // nothing to do
	
	bm->indicators = NULL;
	bm->unscaled_rates = NULL;
	bm->map = NULL;
	
	
	return bm;
}

// consider the value is within bounds!
void _set_strickclock( BranchModel *bm, const int index, const double value ){
	Parameters_set_value(bm->rates, 0, value);
}

double _get_strict_clock( BranchModel *bm , Node *node ){
	return Parameters_value(bm->rates, 0);
}

#pragma mark -
// MARK: Local clock

/***************************** LOCAL CLOCK *****************************/

static double _get_local_clock( BranchModel *bm , Node *node );
static void _set_parameter_localclock( BranchModel *bm, const int index, const double value );

/** 
 * Create a LocalClock structure
 * @param tree: a tree
 * @param n: number of local clocks
 */
BranchModel * new_LocalClock( Tree *tree, const int n ){
	Parameters *rates = new_Parameters(n+1);
	char name[50] = "rate";
	for (int i = 0; i <= n; i++) {
		sprintf(name+4, "%d", i);
		Parameters_move(rates, new_Parameter_with_postfix(name, LOCAL_POSTFIX, 0.001, new_Constraint(BRANCHMODEL_RATE_MIN, BRANCHMODEL_RATE_MAX)));
	}
	return new_LocalClock_with_parameters( tree, rates );
}

BranchModel * new_LocalClock_with_parameters( Tree *tree, const Parameters *rates ){
	BranchModel *bm = BranchModel_init(tree, CLOCK_LOCAL);
		
	bm->get = _get_local_clock;
	bm->set = _set_parameter_localclock;
	
	bm->indicators = bvector(Tree_node_count(tree));
	bm->unscaled_rates = dvector(Tree_node_count(tree));
	
	bm->rates = new_Parameters(Parameters_count(rates));
	Parameters_add_parameters(bm->rates, rates);
	
	for (int i = 0; i < Tree_node_count(tree); i++) {
		bm->indicators[i]     = false;
		bm->unscaled_rates[i] = 0.0;
	}
	
	localclock_set_random_clock_indicators( bm, Parameters_count(rates)-1);
	
	bm->map = new_DiscreteParameter_with_postfix("map", LOCAL_POSTFIX, Tree_node_count(tree) );
	localclock_rebuild_map( bm );
	
	bm->need_update = false;
	
	return bm;
}

// If rates are not present in the tree node descriptions their values are NAN
BranchModel * new_LocalClock_from_tree( Tree *tree ){
    
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	int nNodes = Tree_node_count(tree);	
	
	bool *indicators = bvector( nNodes );
	Vector *vrates = new_Vector(2);
	iVector *pre_order = new_iVector(2);
	double rate = 0;
	
    // The main global rate from the root
	if ( String_contains_str( Tree_root(tree)->right->info, "local=1") ) {
		rate = Node_get_double_from_info(Tree_root(tree)->left, "rate=");
	}
	else {				
		rate = Node_get_double_from_info(Tree_root(tree)->right, "rate=");
	}
	
	Vector_push(vrates, rate);			
	iVector_push( pre_order, 0);	
	
	for ( int i = 0; i < nNodes-1; i++ ) {
		if ( String_contains_str( nodes[i]->info, "local=1") ) {
			indicators[ Node_id(nodes[i]) ] = true;
			rate = Node_get_double_from_info(nodes[i], "rate=");
			Vector_push(vrates, rate);
			iVector_push( pre_order, Node_id(nodes[i]));
		}
	}
	
	BranchModel *bm = new_LocalClock(tree, Vector_length(vrates)-1);
	
	localclock_set_indicators(bm, indicators);
	
	// reorder rates because indexes in map are set in preorder and rates were added in postorder
	Vector_sort_from_iVector(vrates, pre_order);
	
	for ( int i = 0; i < Vector_length(vrates); i++ ) {
		Parameters_set_value(bm->rates, i, Vector_at(vrates, i));
	}
    
    /*for ( int i=0; i < nNodes; i++ ) {
        printf("%s %d %d %f\n", nodes[i]->name, bm->map[nodes[i]->id], indicators[nodes[i]->id], bm->get(bm, nodes[i]));
    }*/
	
	free(indicators);
	free_Vector(vrates);
	free_iVector(pre_order);
	
	return bm;
}


void LocalClock_set_number_of_clocks( BranchModel *bm, const int nLocalClocks ){
	int count = Parameters_count(bm->rates)-1;
	if ( count == nLocalClocks ) return;
	
	if ( count > nLocalClocks ) {
		for (int i = 0; i < count - nLocalClocks; i++) {
			Parameters_pop( bm->rates );
		}
	}
	else{
		char name[50] = "rate";
		for (int i = count; i < nLocalClocks; i++) {
			sprintf(name+4, "%d", (i+1) );
			double rate = Parameters_value(bm->rates, count);
			Parameter *p = new_Parameter_with_postfix(name, LOCAL_POSTFIX, rate, Parameters_constraint(bm->rates, 0)); // getting the common constraint from Parameters not Parameter
			Parameters_move(bm->rates, p);
		}
	}
	
	for (int i = 1; i < Parameters_count(bm->rates); i++) {
		Parameters_set_value(bm->rates, i, Parameters_value(bm->rates, 0) );
	}
	
	localclock_set_random_clock_indicators(bm, nLocalClocks);
	
	localclock_rebuild_map(bm);
}

void localclock_set_random_clock_indicators( BranchModel *bm, const int nLocalClocks ){
	assert( nLocalClocks < Tree_node_count(bm->tree) );
	memset(bm->indicators, false, Tree_node_count(bm->tree) * sizeof(bool) );
	int count = 0;
    int root_id = Node_id(Tree_root(bm->tree));
    
	while ( count != nLocalClocks ) {
		int b = random_int( Tree_node_count(bm->tree)-1 );
        // we never assign a local clock to the root
		if ( root_id != b && !bm->indicators[b] ) {
			count++;
			bm->indicators[b] = true;
		}
	}
	
}

// positions contains the node ids for local clocks
void localclock_set_indicators2( BranchModel *bm, const unsigned int *positions ){
	memset(bm->indicators, false, Tree_node_count(bm->tree) * sizeof(bool) );
	for (int i = 0; i < Parameters_count(bm->rates)-1; i++) {
		bm->indicators[ positions[i] ] = true;
	}
	localclock_rebuild_map( bm );
}


void localclock_set_indicators( BranchModel *bm, const bool *indicators ){
	memcpy( bm->indicators, indicators, Tree_node_count(bm->tree) * sizeof(bool) );
	localclock_rebuild_map( bm );
}

// The map is indexed by node id. Each element contains the index to a rate.
// Using preorder so a node inherits the rate from its parent.
// In postorder we would have to do another pass to propagate the index down the tree
// If indicator == 1 then then the node starts a local clock <-> inheritable rate (rate switch coincide with speciation)
void LocalClock_indicator_to_map( const bool *indicators, unsigned int *map, Node *node, int *index ){
	if( node == NULL ) return;
	
	if( indicators[ Node_id(node) ] ){
		(*index)++;
		map[ Node_id(node) ] = *index;
	}
	else if( Node_isroot(node) ){
		map[ Node_id(node) ] = 0; // this is the root, should not change, SHOULD check when I create, mutate and mate individuals
	}
	else{
		map[ Node_id(node) ] = map[ Node_id(Node_parent(node)) ];
	}
	
	
	LocalClock_indicator_to_map(indicators, map, Node_left(node),  index );
	LocalClock_indicator_to_map(indicators, map, Node_right(node), index );
}

// The map is indexed by node id. Elements contain the index of the rate.
// Using preorder so a node inherits the rate from its parent.
// In postorder we would have to do another pass to propagate the index down the tree
// If indicator == 1 then then its children node start the same local clock <-> rate switch occur just before speciation (branching)
// For external nodes it is different. If indicator is 1 it has its own rate, otherwise 2 sister taxa would always have the same rate.
void LocalClock_indicator_to_map_not_inheritable( const bool *indicators, unsigned int *map, Node *node, int *index ){
	if( node == NULL ) return;
	
    if( Node_isroot(node) ){
		map[ Node_id(node) ] = 0; // this is the root, should not change, SHOULD check when I create, mutate and mate individuals
	}
    // Extremely hackish
    else if( Node_isleaf(node) && indicators[ Node_id(node) ] ){
        (*index)++;
        map[ Node_id(node) ] = *index;
    }
	else if( indicators[ Node_id(Node_parent(node)) ] ){
		if(Node_parent(node)->left == node ){
            (*index)++;
            map[ Node_id(node) ] = *index;
        }
        else {
            map[ Node_id(node) ] = map[ Node_id(node->parent->left) ];
        }
	}
	else{
		map[ Node_id(node) ] = map[ Node_id(Node_parent(node)) ];
	}
	
	
	LocalClock_indicator_to_map_not_inheritable(indicators, map, node->left,  index );
	LocalClock_indicator_to_map_not_inheritable(indicators, map, node->right, index );
}

void localclock_rebuild_map( BranchModel *bm ){
	int index = 0;
	LocalClock_indicator_to_map( bm->indicators, bm->map->values, Tree_root(bm->tree), &index);
}

/**
 * @param index id of the node
 */
double _get_local_clock( BranchModel *bm , Node *node ){
	//if( index == Tree_node_count(bm->tree)-1 ) error("_get_local_clock: should not be called on the root\n");
	return Parameters_value( bm->rates, bm->map->values[ Node_id(node) ] );
}

/**
 * @param index index of the rate in the list (not the index of the branch!)
 * @param value value of the parameter
 */
void _set_parameter_localclock( BranchModel *bm, const int index, const double value ){
	Parameters_set_value(bm->rates, index, value);
}

// indexes contains the node ids where there are local clocks
void LocalClock_get_indexes( const BranchModel *bm,  unsigned *indexes){
    Node **nodes = Tree_nodes(bm->tree);
	int count = 0;
	for (int i = 0; i < Tree_node_count(bm->tree); i++) {
		if ( bm->indicators[i] ) {
			indexes[count++] = Node_id(nodes[i]);
		}
	}
}

#pragma mark -
// MARK: DiscreteClock

/***************************** DISCRETE CLOCK *****************************/

static double _get_DiscreteClock2( BranchModel *bm , Node *node );
static double _get_DiscreteClock( BranchModel *bm , Node *node );
static void _set_DiscreteClock( BranchModel *bm, const int index, const double value );

BranchModel * new_DiscreteClock2( Tree *tree, const int n ){
	Parameters *rates = new_Parameters(n);
	char name[50] = "rate";
	for (int i = 0; i < n; i++) {
		sprintf(name+4, "%d", i);
		Parameters_move(rates, new_Parameter_with_postfix(name, DISCRETE_POSTFIX, .01, new_Constraint(BRANCHMODEL_RATE_MIN, 100)));
	}
    BranchModel *bm = new_DiscreteClock_with_parameters( tree, rates );
    bm->get = _get_DiscreteClock2;
	return bm;
}

BranchModel * new_DiscreteClock( Tree *tree, const int n ){
	Parameters *rates = new_Parameters(n);
	char name[50] = "rate";
	for (int i = 0; i < n; i++) {
		sprintf(name+4, "%d", i);
		Parameters_move(rates, new_Parameter_with_postfix(name, DISCRETE_POSTFIX, .01, new_Constraint(BRANCHMODEL_RATE_MIN, 100)));
	}
	return new_DiscreteClock_with_parameters( tree, rates );
}

BranchModel * new_DiscreteClock_with_parameters( Tree *tree, const Parameters *rates ){
	BranchModel *bm = BranchModel_init(tree, CLOCK_DISCRETE);
		
	bm->get = _get_DiscreteClock;
	bm->set = _set_DiscreteClock;
	
	bm->rates = new_Parameters(Parameters_count(rates));
	Parameters_add_parameters(bm->rates, rates);
	
	bm->map = new_DiscreteParameter_with_postfix("map", DISCRETE_POSTFIX, Tree_node_count(tree) );
	
	DiscreteClock_set_random_branch_assigment( bm );
	
	bm->indicators = NULL;
	bm->unscaled_rates = NULL;
	bm->need_update = false;
	
	return bm;
}


BranchModel * new_DiscreteClock_from_tree( Tree *tree ){
    
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	int nNodes = Tree_node_count(tree);
	unsigned *map = uivector(nNodes);
	int nclasses = 0;
    
    
	for ( int i = 0; i < nNodes-1; i++ ) {
		if ( String_contains_str( nodes[i]->info, "class=") ) {
			int class = Node_get_int_from_info(nodes[i], "class=");
			nclasses = imax(class, nclasses);
			map[ Node_id(nodes[i]) ] = class;
		}
        else {
            error("Missing class in new_DiscreteClock_from_tree\n");
        }
	}
    
    int root_id = Node_id(Tree_root(tree));
    if ( root_id == 0 ) {
        map[0] = map[1];
    }
    else {
        map[root_id] = map[root_id-1];
    }
    
	nclasses++;
	
	double *rates = dvector(nclasses);
	
	for ( int i = 0; i < nNodes-1; i++ ) {
		if ( String_contains_str( nodes[i]->info, "class=") ) {
			int class   = Node_get_int_from_info(nodes[i], "class=");
			rates[class] = Node_get_double_from_info(nodes[i], "rate=");
		}
	}
	
	// at this point the indexes can be unordered in postorder(e.g. 1203 instead of 0123)
    // I can't remember why I should care about the ordering at this point. It is important in the GA
    // but it should be taken care of inside the GA.
	int *v =  ivector(nclasses);
	
	uivector_canonical(map, nNodes-1, v, nclasses);
	
	dvector_sort_from_ivector(rates, v, nclasses);
	free(v);
	
	BranchModel *bm = new_DiscreteClock(tree, nclasses);
	
	DiscreteClock_set_classes(bm, map);
	BranchModel_vector_to_rates(bm, rates);
	free(map);
	free(rates);
	
	return bm;
}

BranchModel * new_DiscreteClock_from_LocalClock_tree( Tree *tree ){
	//TODO: indicators and map are indexed by id not postorder_idx
    fprintf(stderr, "Replace postorder_idx with id %s %d\n",__FILE__,__LINE__);
    exit(1);
	Node **temp_nodes = Tree_get_nodes(tree, POSTORDER);
	int nNodes = Tree_node_count(tree);	
	
	unsigned *map = uivector(nNodes);
	
	bool *indicators = bvector( nNodes );
	Vector *vrates = new_Vector(2);
	iVector *pre_order = new_iVector(2);
	double rate = 0;
	
	if ( String_contains_str( Tree_root(tree)->right->info, "local=1") ) {
		rate = Node_get_double_from_info(Tree_root(tree)->left, "rate=");
	}
	else {				
		rate = Node_get_double_from_info(Tree_root(tree)->right, "rate=");
	}
	
	Vector_push(vrates, rate);			
	iVector_push( pre_order, 0);	
	
	for ( int i = 0; i < nNodes-1; i++ ) {
		if ( String_contains_str( temp_nodes[i]->info, "local=1") ) {
			indicators[ temp_nodes[i]->postorder_idx ] = true;
			rate = Node_get_double_from_info(temp_nodes[i], "rate=");
			Vector_push(vrates, rate);
			iVector_push( pre_order, temp_nodes[i]->preorder_idx);
		}
	}
	
	int index = 0;
	
	// create map from indicators
	LocalClock_indicator_to_map(indicators, map, Tree_root(tree), &index);
	
	// reorder rates because indexes are set in preorder and rates were added in postorder
	Vector_sort_from_iVector(vrates, pre_order);
	
	// at this point the indexes can be unordered in postorder(e.g. 1203 instead of 0123)
	int nclasses = Vector_length(vrates);
		
	int *v =  ivector(nclasses);
	uivector_canonical(map, nNodes-1, v, nclasses);
	
	Vector_sort_from_ivector(vrates, v);
	free(v);
	
	BranchModel *bm_discrete = new_DiscreteClock(tree, nclasses);
	
	DiscreteClock_set_classes(bm_discrete, map);
	
	for ( int i = 0; i < nclasses; i++ ) {
		bm_discrete->set( bm_discrete, i, Vector_at(vrates, i) );
	}
	
	
	free_iVector(pre_order);
	free_Vector(vrates);
	free(indicators);
	free(map);

	return bm_discrete;
}

BranchModel * new_DiscreteClock_from_LocalClock( const BranchModel *localBm ){
	
	Parameters *newRates = clone_Parameters(localBm->rates);
	
	BranchModel *bm = new_DiscreteClock_with_parameters( localBm->tree, newRates);
	DiscreteClock_set_classes(bm, localBm->map->values);
	bm->id = localBm->id+1;
	return bm;
}

void DiscreteClock_set_number_of_rate_classes( BranchModel *bm, const int nClasses ){
	int count = Parameters_count(bm->rates);
	if ( count == nClasses ) return;
	if ( count > nClasses ) {
		for (int i = 0; i < count - nClasses; i++) {
			Parameters_pop( bm->rates );
		}
	}
	else{
		char name[50] = "rate";
		double rate = Parameters_value(bm->rates, count-1);
		//fprintf(stderr, "Add rates before %d after %d\n",count,n);
		for (int i = count; i < nClasses; i++) {
			sprintf(name+4, "%d", i );
			Parameters_move(bm->rates, new_Parameter_with_postfix(name, DISCRETE_POSTFIX, rate, Parameters_constraint(bm->rates, -1)));
		}
		//fprintf(stderr, "Real number %d\n", Parameters_count(bm->rates)-1);
	}
	
	/*for (int i = 1; i < Parameters_count(bm->rates); i++) {
		//bm->rates->list[i]->value = 1.;
		Parameters_set_value(bm->rates, i, Parameters_value(bm->rates, 0) );
	}*/
	
	DiscreteClock_set_random_branch_assigment(bm);
	
}


void DiscreteClock_set_random_branch_assigment( BranchModel *bm ){
	for ( int i = 0; i < Tree_node_count(bm->tree); i++) {
		bm->map->values[i] = random_int(Parameters_count(bm->rates)-1);
	}
}

static double calculateScaleFactor( BranchModel *bm ) {
    double sumTime = 0.0;
    double sumDistance = 0.0;
    Node **nodes = Tree_nodes(bm->tree);
    
    for (int i = 0; i < Tree_node_count(bm->tree); i++) {
        if( Node_isroot(nodes[i]) ) continue;
        
        double times = Node_time_elapsed(nodes[i]);
        double distances = times * Parameters_value(bm->rates, bm->map->values[i]);
        
        sumTime     += times;
        sumDistance += distances;
    }
    
    return sumTime / sumDistance;
}

/**
 * @param index id of the node
 */
double _get_DiscreteClock( BranchModel *bm , Node *node ){
    return Parameters_value(bm->rates, bm->map->values[ Node_id(node) ] );
}

/**
 * @param index id of the node
 */
double _get_DiscreteClock2( BranchModel *bm , Node *node ){
    bm->scalefactor = calculateScaleFactor(bm);
    
    //printf("scaleFactor %f\n", scaleFactor);
	return Parameters_value(bm->rates, bm->map->values[ Node_id(node) ] )*bm->scalefactor*Parameters_value(bm->rates, Parameters_count(bm->rates)-1);
}

void _set_DiscreteClock( BranchModel *bm, const int index, const double value ){
	//fprintf(stderr, "\n%s: _set_localclock: %s %f %f %d\n", __FILE__, bm->rates->list[index]->name, bm->rates->list[index]->value, value, index);
	Parameters_set_value(bm->rates, index, value);
}


void DiscreteClock_set_classes( BranchModel *bm, const unsigned int *classes ){
	memcpy(bm->map, classes, Tree_node_count(bm->tree) * sizeof(unsigned int) );
}

#pragma mark -
// MARK: Relaxed Clock

/***************************** RELAXED CLOCK *****************************/

static double _get_RelaxedClock( BranchModel *bm , Node *node );
static void _set_RelaxedClock( BranchModel *bm, const int index, const double value );

static Parameters * _lognormal_create_parameters( const double logmean, const double logsigma );
static Parameters * _exponetial_create_parameters( const double lambda );
static Parameters * _log_spaced_create_parameters( const double center );

BranchModel * new_RelaxedClock( Tree *tree, const relaxed_clock type, const int n, ... ){
	va_list ap;
	
	va_start(ap,n);
	Parameters *ps = NULL;
	
	switch ( type ) {
		case RELAXED_LOGNORMAL:{
			double logmean  = va_arg(ap, double);
			double logsigma = va_arg(ap, double);
			ps = _lognormal_create_parameters( logmean, logsigma );
			break;
		}
		case RELAXED_EXPONENTIAL:{
			double lambda = va_arg(ap, double);
			ps = _exponetial_create_parameters(lambda);
			break;
		}
        case RELAXED_DISCRETE:{
			double center  = va_arg(ap, double);
            ps = _log_spaced_create_parameters(center);
            break;
        }
		default:
			break;
	}
	
	va_end(ap);
	return new_RelaxedClock_with_parameters( tree, ps, type );
}

BranchModel * new_RelaxedClock_with_parameters( Tree *tree, const Parameters *params, const relaxed_clock type ){
	BranchModel *bm = BranchModel_init(tree, CLOCK_RELAXED);
	
	bm->type = type;
	
	bm->get = _get_RelaxedClock;
	bm->set = _set_RelaxedClock;
	
	bm->rates = new_Parameters(Parameters_count(params));
	Parameters_add_parameters(bm->rates, params);
	
	bm->map = new_DiscreteParameter_with_postfix("map", RELAXED_POSTFIX, Tree_node_count(tree) );
	bm->unscaled_rates = dvector(Tree_node_count(bm->tree));
	
	
	RelaxedClock_set_random_branch_assigment( bm );
	
	bm->indicators = NULL;
	bm->need_update = true;
	
	return bm;
}

BranchModel * new_RelaxedClock_from_tree( Tree *tree, double center ){
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	int nNodes = Tree_node_count(tree);
    
    //FIXME: rate parameter
	BranchModel *bm = new_RelaxedClock(tree, RELAXED_DISCRETE, 1, center);
	for ( int i = 0; i < nNodes-1; i++ ) {
		bm->map->values[i] = i;
		bm->unscaled_rates[i] = Node_get_double_from_info(nodes[i], "rate=");
	}
	bm->need_update = false;
	return bm;
}


Parameters * _lognormal_create_parameters( const double logmean, const double logsigma ){
	Parameters *ps = new_Parameters(2);
	Parameters_move(ps, new_Parameter_with_postfix("logmean", RELAXED_POSTFIX, logmean, new_Constraint(-INFINITY, INFINITY) ) );
	Parameters_move(ps, new_Parameter_with_postfix("logsigma", RELAXED_POSTFIX, logsigma, new_Constraint(TINY, INFINITY) ) );
	return ps;
}

Parameters * _exponetial_create_parameters( const double lambda ){
	Parameters *ps = new_Parameters(2);
	Parameters_move(ps, new_Parameter_with_postfix("lambda", RELAXED_POSTFIX, lambda, new_Constraint(TINY, INFINITY) ) );
	return ps;
}

Parameters * _log_spaced_create_parameters( const double center ){
	Parameters *ps = new_Parameters(1);
	Parameters_move(ps, new_Parameter_with_postfix("center", RELAXED_POSTFIX, center, new_Constraint(BRANCHMODEL_RATE_MIN, BRANCHMODEL_RATE_MAX) ) );
	return ps;
}

// set random location but with bm->tree->nNodes == n categories
void RelaxedClock_set_random_branch_assigment( BranchModel *bm ){
	RelaxedClock_set_random_branch_assigment2( bm, Tree_node_count(bm->tree) );
}

// set random location but with cat_count cateogries
void RelaxedClock_set_random_branch_assigment2( BranchModel *bm, const int cat_count ){
	int cat = 0;
	for ( int i = 0 ; i < Tree_node_count(bm->tree); i++,cat++ ) {
		if ( cat == cat_count ) {
			cat = 0;
		}
		bm->map->values[i] = cat;
	}
	
	randomize_uivector( bm->map->values, Tree_node_count(bm->tree) );
}

// array with 0 at the first element followed by exponentialy spaced, followed by log spaced
void _exp_zero_log(BranchModel *bm ){
    int n_lower = Tree_node_count(bm->tree)/2;
    int n_upper = Tree_node_count(bm->tree)/2;
    
    if( !(Tree_node_count(bm->tree) & 1) ){
        n_lower--;
    }
    double magnitude = 10;
    double lower = Parameters_value(bm->rates, 0) / magnitude;
    double upper = Parameters_value(bm->rates, 0) * magnitude;
    
    bm->unscaled_rates[0] = 0.0;
    
    exp_spaced_spaced_vector2(bm->unscaled_rates+1, lower, Parameters_value(bm->rates, 0), n_lower);
    log_spaced_spaced_vector2(bm->unscaled_rates+1+n_lower, Parameters_value(bm->rates, 0), upper, n_upper);
}

void _relaxedclock_calculate_rates( BranchModel *bm ){
	switch ( bm->type) {
		case RELAXED_LOGNORMAL:{
			lognormal_discretize( Parameters_value(bm->rates, 0), Parameters_value(bm->rates, 1), bm->unscaled_rates, Tree_node_count(bm->tree) );
			break;
        }
		case RELAXED_EXPONENTIAL:{
			exponential_discretize( Parameters_value(bm->rates, 0), bm->unscaled_rates, Tree_node_count(bm->tree) );
			break;
        }
        case RELAXED_DISCRETE:{
            int n_lower = Tree_node_count(bm->tree)/2;
            int n_upper = Tree_node_count(bm->tree)-n_lower;
            
            double magnitude = 10;
            double lower = Parameters_value(bm->rates, 0) / magnitude;
            double upper = Parameters_value(bm->rates, 0) * magnitude;

            exp_spaced_spaced_vector2(bm->unscaled_rates, lower, Parameters_value(bm->rates, 0), n_lower);
            log_spaced_spaced_vector2(bm->unscaled_rates+n_lower, Parameters_value(bm->rates, 0), upper, n_upper);
            break;
        }
		default:
			assert(0);
	}
}

double _get_RelaxedClock( BranchModel *bm , Node *node ){
	if( bm->need_update ){
		_relaxedclock_calculate_rates(bm);
		bm->need_update = false;
		
	}
	return bm->unscaled_rates[ bm->map->values[ Node_id(node) ] ];
}

// index is the index of the parameter of the relaxed clock (not the the postorder index)
void _set_RelaxedClock( BranchModel *bm, const int index, const double value ){
    Parameters_set_value(bm->rates, index, value);
    bm->need_update = true;
}


void RelaxedClock_set_classes( BranchModel *bm, const unsigned int *classes ){
	memcpy(bm->map, classes, Tree_node_count(bm->tree) * sizeof(unsigned int) );
}


#pragma mark -
#pragma mark Misc


static void _get_distance_aux( BranchModel *bm, Node *n  ){	
	if( Node_parent(n) != NULL ){
		Node_set_distance( n, bm->get(bm, n) * (Node_height(Node_parent(n)) - Node_height(n) ) );
		if( Node_distance(n) < 0 )
			fprintf(stderr, "%s branch length = %f rate = %f parent height [%s]= %f node height = %f\n", n->name, Node_distance(n), bm->get( bm, n), n->parent->name, Node_height(Node_parent(n)), Node_height(n));
	
	}
	
	if( !Node_isleaf(n) ){
		_get_distance_aux( bm, Node_left(n) );
		_get_distance_aux( bm, Node_right(n) );
	}
}

void infer_distance_from_rate_height( BranchModel *bm ){
	_get_distance_aux( bm, Tree_root(bm->tree) );
}

void print_rate_map( BranchModel *bm ){
	Node **nodes = Tree_get_nodes(bm->tree, POSTORDER);
	fprintf(stderr, "\n--------------------------\n");
	for (int i = 0; i < Tree_node_count(bm->tree)-1; i++) {
		fprintf(stderr, "%s map:%d ", nodes[i]->name, bm->map->values[i]);
		if(bm->indicators != NULL ) fprintf(stderr, "indicator:%d value:%f (param:%f; unscaled:%f)\n", bm->indicators[i], bm->get(bm, nodes[i]), Parameters_value(bm->rates, bm->map->values[ Node_id(nodes[i]) ]), bm->unscaled_rates[i]);
		else fprintf(stderr, "value:%f param:%f\n", bm->get(bm, nodes[i]), Parameters_value(bm->rates, bm->map->values[ Node_id(nodes[i])] ));
	}
	//fprintf(stderr, "%s map:%d indicator:%d value:- (param:%f; unscaled:%f)\n", nodes[bm->tree->nNodes-1]->name, bm->map[bm->tree->nNodes-1], bm->indicators[bm->tree->nNodes-1], bm->rates->list[bm->map[nodes[bm->tree->nNodes-1]->postorder_idx]]->value, bm->unscaled_rates[bm->tree->nNodes-1]);
	fprintf(stderr, "\n==========================\n");
}



void print_compare_bms( const BranchModel *bm1, const BranchModel *bm2){
	printf(" ID %d %d\nNAME %d %d\nNEED UPDATE %d %d\nSCALEFACTOR %f %f\n",bm1->id, bm2->id,bm1->name, bm2->name,bm1->need_update,bm2->need_update,bm1->scalefactor,bm2->scalefactor);
	Parameters_print(bm1->rates);
	Parameters_print(bm2->rates);
	
	print_ivector( (int*)bm1->indicators, Tree_node_count(bm1->tree));
	print_ivector( (int*)bm2->indicators, Tree_node_count(bm2->tree));
	
	print_uivector( bm1->map->values, Tree_node_count(bm1->tree));
	print_uivector( bm2->map->values, Tree_node_count(bm2->tree));
	
	print_dvector( bm1->unscaled_rates, Tree_node_count(bm1->tree));
	print_dvector( bm2->unscaled_rates, Tree_node_count(bm2->tree));
}

void compare_branchmodel( const BranchModel *bm1, const BranchModel *bm2 ){
	assert( bm1->name == bm2->name );
	assert( bm1->id == bm2->id );
	compare_tree( bm1->tree, bm2->tree );
	assert( bm1->need_update == bm2->need_update );
	assert( bm1->scalefactor == bm2->scalefactor );
	compare_parameters(bm1->rates, bm2->rates);
	assert( memcmp( bm1->indicators, bm2->indicators, Tree_node_count(bm1->tree) * sizeof(bool) ) == 0 );
	assert( memcmp( bm1->map, bm2->map, Tree_node_count(bm1->tree) * sizeof(int) ) == 0 );
	assert( memcmp( bm1->unscaled_rates, bm2->unscaled_rates, Tree_node_count(bm1->tree) * sizeof(double) ) == 0 );
}

double BranchModel_mean_rate_scaled( BranchModel *bm ){
    if( bm->name == CLOCK_STRICT) return Parameters_value(bm->rates, 0);
    
	Node **nodes = Tree_nodes(bm->tree);
	double d = 0.0;
	double t = 0.0;
	for ( int i = 0; i < Tree_node_count(bm->tree); i++ ) {
        if( Node_isroot(nodes[i]) ) continue;
		double bt = Node_height( Node_parent(nodes[i]) ) - Node_height( nodes[i] );
		d += bm->get(bm, nodes[i]) * bt;
		t += bt;
	}
	return d/t;
}

double BranchModel_mean_rate_tips_scaled( BranchModel *bm ){
	Node **nodes = Tree_nodes(bm->tree);
	double d = 0.0;
	double t = 0.0;
	for ( int i = 0; i < Tree_node_count(bm->tree); i++ ) {
        if( Node_isroot(nodes[i]) ) continue;
        
        if ( Node_isleaf(nodes[i]) ) {
            double bt = Node_height( Node_parent(nodes[i]) ) - Node_height( nodes[i] );
            d += bm->get(bm, nodes[i]) * bt;
            t += bt;
        }
	}
	return d/t;
}

double BranchModel_mean_rate_internal_scaled( BranchModel *bm ){
	Node **nodes = Tree_nodes(bm->tree);
	double d = 0.0;
	double t = 0.0;
	for ( int i = 0; i < Tree_node_count(bm->tree); i++ ) {
        if( Node_isroot(nodes[i]) ) continue;
        
        if ( !Node_isleaf(nodes[i]) ) {
            double bt = Node_height( Node_parent(nodes[i]) ) - Node_height( nodes[i] );
            d += bm->get(bm, nodes[i]) * bt;
            t += bt;
        }
	}
	return d/t;
}

double BranchModel_mean_rate( BranchModel *bm, double *min, double *max ){
	*min = 100;
	double meanRate = 0;
	for ( int i = 0; i < Parameters_count(bm->rates); i++ ) {
		meanRate += Parameters_value(bm->rates, i);
		*max = dmax(*max, Parameters_value(bm->rates, i));
		*min = dmin(*min, Parameters_value(bm->rates, i));
	}
	return meanRate/Parameters_count(bm->rates);
}

// return rate correlation between adjacent branches
double BranchModel_correlation( BranchModel *bm ) {
	Tree *tree = bm->tree;
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	double *x = dvector(Tree_node_count(tree)-3);
	double *y = dvector(Tree_node_count(tree)-3);
	int index = 0;
	for (int i = 0; i < Tree_node_count(tree)-2; i++) {
		Node *parent = Node_parent(nodes[i]);
		if ( !Node_isroot(parent) ) {
			x[index] = bm->get(bm, nodes[i]);
			y[index] = bm->get(bm, parent);
			index++;
		}
	}
	double cor = correlation(x,y,Tree_node_count(tree)-3);
	free(x);
	free(y);
	return cor;
}

// return the correlation between branch length estimated with and without rate
double BranchModel_correlation_distance( BranchModel *bm ){
	Node **nodes = Tree_nodes(bm->tree);
    double *b = dvector(Tree_node_count(bm->tree)-2);
    double *b2 = dvector(Tree_node_count(bm->tree)-2);
    
    for (int i = 0; i < Tree_node_count(bm->tree); i++) {
        
        if( Node_isroot(nodes[i]) || ( Node_isroot(nodes[i]->parent) && nodes[i]->parent->right == nodes[i] ) ) continue;
        
        b[i]  = Node_distance(nodes[i]);
        b2[i] = bm->get(bm, nodes[i]) * Node_time_elapsed(nodes[i]);
        
        if( Node_isroot( Node_parent(nodes[i]) ) ){
            Node *sibling = Node_sibling(nodes[i]);
            b[i]  += Node_distance(sibling);
            b2[i] += bm->get(bm, sibling) * Node_time_elapsed(sibling);
        }
        //printf("%f %f\n", b[i],b2[i]);
    }
    double cor = correlation(b, b2,Tree_node_count(bm->tree)-2);
    
    //fprintf(stderr, "Correlation = %f  (%s)\n\n", cor, argv[2]);
    
    free(b);
    free(b2);
    return cor;
}

static void grab_distances( BranchModel *bm, Node *n, double *x, double *y, int *index  ){
	if( !Node_isroot(n) ){
		double bl1 = Node_distance(n);
		if(bl1 < 0 ){
			fprintf(stderr, "grab_distances: %s branch length = %f", n->name, bl1);
		}
		
		double bl2 = bm->get(bm, n) * (Node_height(Node_parent(n)) - Node_height(n) );
		if(bl2 < 0 )
			fprintf(stderr, "grab_distances: %s branch length = %f rate = %f height = %f - parent height [%s]= %f\n", n->name, bl2, bm->get(bm, n), Node_height(n), n->parent->name, Node_height(Node_parent(n)));
		x[*index] = bl1;
		y[*index] = bl2;
		(*index)++;
	}
	
	if( !Node_isleaf(n) ){
		grab_distances( bm, Node_left(n),x, y, index );
		grab_distances( bm, Node_right(n), x, y,index );
	}
}

void BranchModel_check_outliers( BranchModel * bm ){
	Tree *tree = bm->tree;
	double *x = dvector(Tree_node_count(tree)-2);
	double *y = dvector(Tree_node_count(tree)-2);
	int index = 0;
	grab_distances(bm, Tree_root(tree), x, y, &index);
	fprintf(stdout, "rate-free,rate\n");
	for (int i = 0; i < Tree_node_count(tree)-2; i++) {
		fprintf(stdout, "%f,%f\n", x[i],y[i]);
	}
	free(x);
	free(y);
}

void BranchModel_to_distance( BranchModel *bm ){
	Tree *tree = bm->tree;
	Node **nodes = Tree_nodes(tree);
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
        if( Node_isroot(nodes[i]) ) continue;
		Node_set_distance(nodes[i], (Node_height(nodes[i]->parent)-Node_height(nodes[i])) );
	}
}

