// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <math.h>

#include "treestat.h"

#include "statistics.h"
#include "tree.h"
#include "node.h"
#include "matrix.h"
#include "mstring.h"


double TreeStat_rate_correlation( const Tree *tree ) {
	Node **nodes = Tree_get_nodes( (Tree *)tree, POSTORDER);
	double *x = dvector(Tree_node_count(tree)-2);
	double *y = dvector(Tree_node_count(tree)-2);
	int index = 0;
	for (int i = 0; i < Tree_node_count(tree)-2; i++) {
		Node *parent = Node_parent(nodes[i]);
		if ( !Node_isroot(parent) ) {
			x[index] = Node_get_double_from_info( nodes[i], "rate=");
			y[index] = Node_get_double_from_info( nodes[i]->parent, "rate=");
			index++;
		}
	}
	double cor = correlation(x,y,Tree_node_count(tree)-2);
	free(x);
	free(y);
	return cor;
}

double TreeStat_mean_rate_scaled( const Tree *tree ){
	Node **nodes = Tree_get_nodes((Tree *)tree, POSTORDER);
	double d = 0.0;
	double t = 0.0;
	for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
		double bt = Node_height( Node_parent(nodes[i]) ) - Node_height( nodes[i] );
		d += Node_get_double_from_info( nodes[i], "rate=") * bt;
		t += bt;
	}
	return d/t;
}

double TreeStat_mean_rate_tips_scaled( const Tree *tree ){
	Node **nodes = Tree_get_nodes((Tree *)tree, POSTORDER);
	double d = 0.0;
	double t = 0.0;
	for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
        if( Node_isleaf(nodes[i]) ){
            double bt = Node_height( Node_parent(nodes[i]) ) - Node_height( nodes[i] );
            d += Node_get_double_from_info( nodes[i], "rate=") * bt;
            t += bt;
        }
	}
	return d/t;
}

double TreeStat_mean_rate_internal_scaled( const Tree *tree ){
	Node **nodes = Tree_get_nodes((Tree *)tree, POSTORDER);
	double d = 0.0;
	double t = 0.0;
	for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
        if( !Node_isleaf(nodes[i]) ){
            double bt = Node_height( Node_parent(nodes[i]) ) - Node_height( nodes[i] );
            d += Node_get_double_from_info( nodes[i], "rate=") * bt;
            t += bt;
        }
	}
	return d/t;
}

double TreeStat_mean_rate( const Tree *tree, double *min, double *max ){
	Node **nodes = Tree_get_nodes((Tree *)tree, POSTORDER);
	*min = INFINITY;
	double meanRate = 0;
	for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
		meanRate += Node_get_double_from_info( nodes[i], "rate=" );
		*max = fmax( *max, Node_get_double_from_info( nodes[i], "rate=") );
		*min = fmin( *min, Node_get_double_from_info( nodes[i], "rate=") );
	}
	return meanRate/Tree_node_count(tree)-1;
}
