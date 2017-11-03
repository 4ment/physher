/*
 *  roottotip.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 8/6/12.
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

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "tree.h"
#include "parameters.h"
#include "optimizer.h"
#include "lm.h"
#include "utils.h"
#include "matrix.h"
#include "statistics.h"

typedef struct BData{
	double *X;
	double *Y;
	Tree *tree;
	double bl_tot;
    
    int n; // number of data points
    int index; // index of the node being optimized
    bool *indicators;
}BData;

// order decreasing correlation
int cmpFunc(const void *p1, const void *p2){
	double **a = (double **)p1;
	double **b = (double **)p2;
	if( a[0][10] < b[0][10]) return 1;
	else if(a[0][10] == b[0][10]) return 0;
	else return -1;
}

// order increasing residuals
int cmpFuncResiduals(const void *p1, const void *p2){
	double **a = (double **)p1;
	double **b = (double **)p2;
	if( a[0][12] > b[0][12]) return 1;
	else if(a[0][12] == b[0][12]) return 0;
	else return -1;
}

void compute_tips_root_distances( double *y, Tree *tree ){
	Node **tips  = get_tips( tree, POSTORDER );
	for ( int k = 0; k < Tree_tip_count(tree); k++ ){
		y[k] = Tree_distance_to_root(tips[k]);
	}
	free(tips);
}

void compute_tips_internal_distances( const Node *root, Node *node, double *Y, int *pos, bool *indicators){
    if( indicators[node->postorder_idx] && node != root ){
        // we could take the age of the estimated internal node
    }
    else {
        if( !Node_isleaf(node) ){
            compute_tips_internal_distances(root, node->left,  Y, pos, indicators);
            compute_tips_internal_distances(root, node->right, Y, pos, indicators);
        }
        else {
            Y[*pos] = Tree_distance_to_root(node);// only up  to root!!
            (*pos)++;
        }
    }
}


void adjust_root( Node *node, double value, double total ){
	Node_set_distance(node->left, value);
	Node_set_distance(node->right, total-value);
}

double objectiveFunction( Parameters *p, double *grad, void *data ){
	BData *d = (BData*)data;
	adjust_root(Tree_root(d->tree), Parameters_value(p, 0), d->bl_tot);
	compute_tips_root_distances(d->Y, d->tree);
	return -correlation(d->X, d->Y, Tree_tip_count(d->tree));
}

double objectiveFunction2( Parameters *p, double *grad, void *data ){
	BData *d = (BData*)data;
	adjust_root(Tree_root(d->tree), Parameters_value(p, 0), d->bl_tot);
	compute_tips_root_distances(d->Y, d->tree);
    double slope,intercept;
    regression( d->X, d->Y, Tree_tip_count(d->tree), &slope, &intercept );
    
    double sumsqr = 0;
    for ( int i = 0; i < Tree_tip_count(d->tree); i++ ) {
        double temp = d->Y[i] - (slope * d->X[i] + intercept);
        sumsqr += temp * temp;
    }
	return sumsqr;
}

// for clustering local clock
double objectiveFunction3( Parameters *p, double *grad, void *data ){
	BData *d = (BData*)data;
    Node *node = Tree_get_node(d->tree, POSTORDER, d->index);
    
	adjust_root(node, Parameters_value(p, 0), d->bl_tot);
    int pos = 0;
	compute_tips_internal_distances(node, node, d->Y, &pos, d->indicators);
    double slope,intercept;
    regression( d->X, d->Y, d->n, &slope, &intercept );
    
    double sumsqr = 0;
    for ( int i = 0; i < d->n; i++ ) {
        double temp = d->Y[i] - (slope * d->X[i] + intercept);
        sumsqr += temp * temp;
    }
	return sumsqr;
}


// If it is going backward it has to start from 0, no restriction when forward as it is converted to backward
double * lm_tree( Tree *tree, bool forward, bool use_correlation ){
    
	Node *nodeRoot = Tree_root(tree);
	double slope, intercept, cor, lowerx, upperx;
	
	Optimizer *opt = new_Optimizer( OPT_BRENT );
	BData *opt_data = (BData*)malloc( sizeof(BData) );
	assert(opt_data);
	opt_data->tree = NULL;
	opt_data->X = dvector(Tree_tip_count(tree));
	opt_data->Y = dvector(Tree_tip_count(tree));
    opt_data->index = 0;
    opt_data->n = Tree_tip_count(tree);
    opt_data->indicators = NULL;
	opt_set_data(opt, opt_data);
    
    if( use_correlation ){
        opt_set_objective_function(opt, objectiveFunction);
    }
    else {
        opt_set_objective_function(opt, objectiveFunction2);
    }
    
	Parameters *ps = new_Parameters(1);
	Parameters_move( ps, new_Parameter("root.pos", 0, new_Constraint(0, 0)) );
	
	Node **tips  = get_tips( tree, POSTORDER );
	
    
    if ( forward ){
        for ( int i = 0; i < Tree_tip_count(tree); i++ ) {
            opt_data->X[i] = tips[i]->time;
        }
    }
    else{
        double max = 0.0;
		for ( int i = 0; i < Tree_tip_count(tree); i++ ) {
			max = dmax(max, tips[i]->time);
		}
        
        for ( int i = 0; i < Tree_tip_count(tree); i++ ) {
            opt_data->X[i] = max-tips[i]->time;
        }
    }
	
	
	opt_data->bl_tot = Node_distance(Tree_root(tree)->left)+Node_distance(Tree_root(tree)->right);
	opt_data->tree = tree;
	Parameters_set_upper(ps, 0, opt_data->bl_tot);
	Parameters_set_value(ps, 0, opt_data->bl_tot/2.0);
	
	
	double fret = 0;
	double status = opt_optimize( opt, ps, &fret);
	if( status == OPT_ERROR ) error("root.pos No SUCCESS!!!!!!!!!!!!\n");
	
	double min = Parameters_value(ps, 0);
    
	compute_tips_root_distances( opt_data->Y, tree );
	
	cor = correlation( opt_data->X, opt_data->Y, Tree_tip_count(tree) );
	regression( opt_data->X, opt_data->Y, Tree_tip_count(tree), &slope, &intercept );
	
	double ci_intercept  = CI_intercept( opt_data->X, opt_data->Y, Tree_tip_count(tree), slope, intercept, 0.975 );
	double ci_slope      = CI_slope( opt_data->X, opt_data->Y, Tree_tip_count(tree), slope, intercept, 0.975 );
	CI_xIntercept( opt_data->X, opt_data->Y, Tree_tip_count(tree), slope, intercept, &lowerx, &upperx,0.975);
	
	double *result    = (double *)malloc(sizeof(double)*13 );
    assert(result);
	
	result[0] = Node_id(nodeRoot);
	
	result[1] = slope;
	result[2] = slope - ci_slope;
	result[3] = slope + ci_slope;
	
	result[4] = intercept;
	result[5] = intercept - ci_intercept;
	result[6] = intercept + ci_intercept;
	
	result[7] = -intercept/ slope;
	result[8] = upperx;
	result[9] = lowerx;
	
	result[10] = cor;
	result[11] = min;
	result[12] = 0.0;//Parameter_value(nodeRoot->distance);
    
    for ( int i = 0; i < Tree_tip_count(tree); i++ ) {
        double temp = opt_data->Y[i] - (slope * opt_data->X[i] + intercept);
        result[12] += temp * temp;
    }
    //result[12] /= Tree_tip_count(tree);
    
	free(tips);
	free_Parameters(ps);
	free(opt_data->X);
	free(opt_data->Y);
	opt_data->tree = NULL;
	free(opt_data);
	free_Optimizer(opt);
	
	return result;
}


double ** max_lm_tree( Tree *tree, bool forward, bool use_correlation ){
    double **results = (double **)malloc( (Tree_node_count(tree)-1) * sizeof(double *) );
    assert(results);
    Node **nodes = Tree_nodes(tree);
    
    int count = 0;
    #pragma omp parallel for
    for( int i = 0; i < Tree_node_count(tree); i++ ){
        if( Node_isroot(nodes[i]) ) continue;
        
        Tree *atree = clone_Tree(tree);

        Tree_reroot(atree, Tree_node(atree, i));
        
        double *result = lm_tree( atree, forward, use_correlation );
#pragma omp critical
{
        results[count]= result;
        results[count][0] = i ;// this is the node id
        count++;
}

        free_Tree(atree);
    }
    
    if( use_correlation ){
        qsort(results, Tree_node_count(tree)-1, sizeof(results[0]), cmpFunc);
    }
    else {
        qsort(results, Tree_node_count(tree)-1, sizeof(results[0]), cmpFuncResiduals);
    }
    
    return results;
}


void _doit( const Node *root, Node *node, double *X, int *pos, bool *indicators){
    if( indicators[node->postorder_idx] && node != root ){
        // we could take the age of the estimated internal node
    }
    else {
        if( !Node_isleaf(node) ){
            _doit(root, node->left,  X, pos, indicators);
            _doit(root, node->right, X, pos, indicators);
        }
        else {
            X[*pos] = node->time;
            (*pos)++;
        }
    }
}



double * lm_tree_cluster( Tree *tree, bool forward, int k ){
    
	Node *nodeRoot = Tree_root(tree);
	double slope, intercept, cor, lowerx, upperx;
    
    int nNodes = Tree_node_count(tree);
	
	Optimizer *opt = new_Optimizer( OPT_BRENT );
	BData *opt_data = (BData*)malloc( sizeof(BData) );
	assert(opt_data);
	opt_data->tree = NULL;
	opt_data->X = dvector(Tree_tip_count(tree));
	opt_data->Y = dvector(Tree_tip_count(tree));
    opt_data->n = 0;
    opt_data->index = nNodes-1;
    opt_data->tree = tree;
    opt_data->indicators = bvector(nNodes);
	opt_set_data(opt, opt_data);
    
    opt_set_objective_function(opt, objectiveFunction3);
    
    
	Parameters *ps = new_Parameters(k);
    StringBuffer * buffer = new_StringBuffer(10);
    for ( int i = 0; i < k; i++ ) {
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "pos%d", i);
        Parameters_move( ps, new_Parameter(buffer->c, 0, new_Constraint(0, 0)) );
    }
    free_StringBuffer(buffer);
	
	
	Node **tips  = get_tips( tree, POSTORDER );
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    
    int *positions = ivector(1);// only for 1 local clock
	
    for ( int i = 0; i < Tree_node_count(tree)-2; i++ ) {
        if( !Node_isleaf(nodes[i]) ){
            Node *node = nodes[i];
            
            if( Node_tip_count(node) > 6 ){
                Tree *tree2 = clone_SubTree(tree, node);
                double *result = lm_tree(tree2, forward, false);
                free_Tree(tree2);
            
                opt_data->n = 0;
                memset(opt_data->indicators, false, nNodes * sizeof(bool) );
                positions[0] = i;
                for (int i = 0; i < k-1; i++) {
                    opt_data->indicators[ positions[i] ] = true;
                }
                
                _doit(Tree_root(tree), Tree_root(tree), opt_data->X, &(opt_data->n), opt_data->indicators);
                printf("in %d out %d tot %d -- ", opt_data->n, Node_tip_count(node), Tree_tip_count(tree));
            
                opt_data->bl_tot = Node_distance(Tree_root(tree)->left)+Node_distance(Tree_root(tree)->right);
                Parameters_set_upper(ps, 0, opt_data->bl_tot);
                Parameters_set_value(ps, 0, opt_data->bl_tot/2.0);
                
                opt_data->index = Tree_root(tree)->postorder_idx;
                
                double fret = 0;
                double status = opt_optimize( opt, ps, &fret);
                if( status == OPT_ERROR ) error("root.pos No SUCCESS!!!!!!!!!!!!\n");
                
                //double min = Parameters_value(ps, 0);
                
                //compute_tips_root_distances( opt_data->Y, tree );
                int pos = 0;
                compute_tips_internal_distances(Tree_root(tree), Tree_root(tree), opt_data->Y, &pos, opt_data->indicators);
                
                double cor = correlation( opt_data->X, opt_data->Y, opt_data->n );
                regression( opt_data->X, opt_data->Y, opt_data->n, &slope, &intercept );
                double sumSqr = 0;
                for ( int i = 0; i < opt_data->n; i++ ) {
                    double temp = opt_data->Y[i] - (slope * opt_data->X[i] + intercept);
                    sumSqr += temp * temp;
                }
                printf("%s slope: %f intercept: %f corr: %f e: %e | slope:%f corr: %f e: %e total e: %e\n", Node_name(node), slope, intercept, cor, sumSqr, result[1], result[10], result[12], (result[12]+sumSqr));
                free(result);
            }
        }
    }
	
    
    if ( forward ){
        for ( int i = 0; i < Tree_tip_count(tree); i++ ) {
            opt_data->X[i] = tips[i]->time;
        }
    }
    else{
        double max = 0.0;
		for ( int i = 0; i < Tree_tip_count(tree); i++ ) {
			max = dmax(max, tips[i]->time);
		}
        
        for ( int i = 0; i < Tree_tip_count(tree); i++ ) {
            opt_data->X[i] = max-tips[i]->time;
        }
    }
	
	
	opt_data->bl_tot = Node_distance(Tree_root(tree)->left)+Node_distance(Tree_root(tree)->right);
	Parameters_set_upper(ps, 0, opt_data->bl_tot);
	Parameters_set_value(ps, 0, opt_data->bl_tot/2.0);
	
	
	double fret = 0;
	double status = opt_optimize( opt, ps, &fret);
	if( status == OPT_ERROR ) error("root.pos No SUCCESS!!!!!!!!!!!!\n");
	
	double min = Parameters_value(ps, 0);
    
	compute_tips_root_distances( opt_data->Y, tree );
	
	cor = correlation( opt_data->X, opt_data->Y, Tree_tip_count(tree) );
	regression( opt_data->X, opt_data->Y, Tree_tip_count(tree), &slope, &intercept );
	
	double ci_intercept  = CI_intercept( opt_data->X, opt_data->Y, Tree_tip_count(tree), slope, intercept, 0.975 );
	double ci_slope      = CI_slope( opt_data->X, opt_data->Y, Tree_tip_count(tree), slope, intercept, 0.975 );
	CI_xIntercept( opt_data->X, opt_data->Y, Tree_tip_count(tree), slope, intercept, &lowerx, &upperx,0.975);
	
	double *result    = (double *)malloc(sizeof(double)*13 );
    assert(result);
	
	result[0] = nodeRoot->postorder_idx;
	
	result[1] = slope;
	result[2] = slope - ci_slope;
	result[3] = slope + ci_slope;
	
	result[4] = intercept;
	result[5] = intercept - ci_intercept;
	result[6] = intercept + ci_intercept;
	
	result[7] = -intercept/ slope;
	result[8] = upperx;
	result[9] = lowerx;
	
	result[10] = cor;
	result[11] = min;
	result[12] = 0.0;//Parameter_value(nodeRoot->distance);
    
    for ( int i = 0; i < Tree_tip_count(tree); i++ ) {
        double temp = opt_data->Y[i] - (slope * opt_data->X[i] + intercept);
        result[12] += temp * temp;
    }
    result[12] /= Tree_tip_count(tree);
    
	free(tips);
	free_Parameters(ps);
	free(opt_data->X);
	free(opt_data->Y);
	opt_data->tree = NULL;
	free(opt_data);
	free_Optimizer(opt);
	
	return result;
}
