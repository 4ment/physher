/*
 *  optimize2.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 8/10/13.
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

#include "tree.h"
#include "node.h"
#include "treelikelihood.h"
#include "optimizer.h"
#include "optimize.h"


typedef struct BrentData2{
	SingleTreeLikelihood **tlk;
    int n;
    function_f f;
	int index_model;
	int index_param;
} BrentData2;



#pragma mark -

static double _standard_loglikelihood_brent( void *data ){
    BrentData2 *mydata = (BrentData2*)data;
    
    double lnl = 0;
    for ( int i = 0; i < mydata->n; i++ ) {
        lnl += mydata->tlk[i]->calculate(mydata->tlk[i]);
    }
    return lnl;
}

// For calculation using lower and upper likelihoods
static double _standard_loglikelihood_upper_brent( void *data ){
    BrentData2 *mydata = (BrentData2*)data;
    
    double lnl = 0;
    for ( int i = 0; i < mydata->n; i++ ) {
        lnl += mydata->tlk[i]->calculate_upper(mydata->tlk[i], Tree_get_node(mydata->tlk[i]->tree, POSTORDER, mydata->index_param));
    }
    return lnl;
}

#pragma mark -

BrentData2 * new_BrentData2( SingleTreeLikelihood **tlk, int n){
	BrentData2 *data = (BrentData2*)malloc( sizeof(BrentData2) );
    assert(data);
	data->index_param = 0;
	data->tlk = tlk;
    data->n = n;
    data->f = _standard_loglikelihood_brent;
	return data;
}

void free_BrentData2( BrentData2 *data ){
	free(data);
	data = NULL;
}

#pragma mark -

static double _brent_optimize_frequency( Parameters *params, double *grad, void *data ){
	BrentData2 *mydata = (BrentData2*)data;
	int index = mydata->index_param;
	SingleTreeLikelihood *stlk = mydata->tlk[mydata->index_model]; // current tlk
	    
    for ( int k = 0; k < mydata->n; k++ ) {
        if( mydata->tlk[k]->sm->m->freqs == stlk->sm->m->freqs ){
            mydata->tlk[k]->sm->m->set_frequency( mydata->tlk[k]->sm->m, Parameters_value(params, 0), index );
            SingleTreeLikelihood_update_all_nodes(mydata->tlk[k]);
        }
    }
	
	//return fabs(stlk->calculate(stlk));
    return fabs(mydata->f(mydata));
}

static double _brent_optimize_relative_rate( Parameters *params, double *grad, void *data ){
	BrentData2 *mydata = (BrentData2*)data;
	int index = mydata->index_param;
	SingleTreeLikelihood *stlk = mydata->tlk[mydata->index_model]; // current tlk
	
    for ( int k = 0; k < mydata->n; k++ ) {
        if( mydata->tlk[k]->sm->m->rates == stlk->sm->m->rates ){
            mydata->tlk[k]->sm->m->set_rate( mydata->tlk[k]->sm->m, Parameters_value(params, 0), index );
            SingleTreeLikelihood_update_all_nodes(mydata->tlk[k]);
        }
    }	
	
	//return fabs(stlk->calculate(stlk));
    return fabs(mydata->f(mydata));
}

static double _brent_optimize_pinv( Parameters *params, double *grad, void *data ){
	BrentData2 *mydata = (BrentData2*)data;
	SingleTreeLikelihood *stlk = mydata->tlk[mydata->index_model]; // current tlk
	
    for ( int k = 0; k < mydata->n; k++ ) {
        if( mydata->tlk[k]->sm->pinv == stlk->sm->pinv ){
            SiteModel_set_pinv( mydata->tlk[k]->sm, Parameters_value(params, 0) ); // does not check
            SingleTreeLikelihood_update_all_nodes(mydata->tlk[k]);
        }
    }
	
	//return fabs(stlk->calculate(stlk));
    return fabs(mydata->f(mydata));
}

static double _brent_optimize_gamma( Parameters *params, double *grad, void *data ){
	BrentData2 *mydata = (BrentData2*)data;
    SingleTreeLikelihood *stlk = mydata->tlk[mydata->index_model]; // current tlk
	
    for ( int k = 0; k < mydata->n; k++ ) {
        if( mydata->tlk[k]->sm->shape == stlk->sm->shape ){
            SiteModel_set_alpha( mydata->tlk[k]->sm, Parameters_value(params, 0) ); // does not check
            SingleTreeLikelihood_update_all_nodes(mydata->tlk[k]);
        }
    }
    
	//return fabs(stlk->calculate(stlk));
    return fabs(mydata->f(mydata));
}

static double _brent_optimize_mu( Parameters *params, double *grad, void *data ){
	BrentData2 *mydata = (BrentData2*)data;
    SingleTreeLikelihood *stlk = mydata->tlk[mydata->index_model]; // current tlk
	
    SiteModel_set_mu( stlk->sm, Parameters_value(params, 0) ); // does not check
    SingleTreeLikelihood_update_all_nodes(stlk);
        
    
	//return fabs(stlk->calculate(stlk));
    return fabs(mydata->f(mydata));
}


static double _brent_optimize_branch_length( Parameters *params, double *grad, void *data ){
	BrentData2 *mydata = (BrentData2*)data;
    
    Node *node = Tree_get_node(mydata->tlk[0]->tree, POSTORDER, mydata->index_param);
    Node_set_distance( node, Parameters_value(params, 0) );
    
    for ( int k = 0; k < mydata->n; k++ ) {
        Node *node = Tree_get_node(mydata->tlk[k]->tree, POSTORDER, mydata->index_param);
        SingleTreeLikelihood_update_one_node( mydata->tlk[k], node);
    }
    return fabs(mydata->f(mydata));
}


#pragma mark -

static double _optimize_brent_relative_rate_all( int tlk_index, Optimizer *opt, BrentData2 *data, Parameters *param ){
	opt_result status = OPT_NEED_CHECK;
    double lnl = 0;
    data->index_model = tlk_index;
    SingleTreeLikelihood *stlk = data->tlk[tlk_index];
    
    for (int i = 0; i < Parameters_count(stlk->sm->m->rates); i++) {
        if( Parameters_fixed( stlk->sm->m->rates, i ) ){
            continue;
        }
        
        data->index_param = i;
        Parameters_set(param, 0, Parameters_at(stlk->sm->m->rates, i));
        
        status = opt_optimize( opt, param, &lnl);
        if( status == OPT_ERROR ) error("OPT.relative_rates No SUCCESS!!!!!!!!!!!!\n");
        
        for ( int k = 0; k < data->n; k++ ) {
            if( data->tlk[k]->sm->m->rates == stlk->sm->m->rates ){
                data->tlk[k]->sm->m->set_rate( data->tlk[k]->sm->m, Parameters_value(param, 0), i );
                SingleTreeLikelihood_update_all_nodes(data->tlk[k]);
            }
        }
    }
	return -lnl;
}

static double _optimize_brent_frequencies_all( int tlk_index, Optimizer *opt, BrentData2 *data, Parameters *param ){
	opt_result status = OPT_NEED_CHECK;
    double lnl = 0;
    data->index_model = tlk_index;
    SingleTreeLikelihood *stlk = data->tlk[tlk_index];
    
    for (int i = 0; i < Parameters_count(stlk->sm->m->freqs); i++) {
        if( Parameters_fixed( stlk->sm->m->freqs, i ) ){
            continue;
        }
        
        data->index_param = i;
        Parameters_set(param, 0, Parameters_at(stlk->sm->m->freqs, i));
        
        status = opt_optimize( opt, param, &lnl);
        if( status == OPT_ERROR ) error("OPT.frequencies No SUCCESS!!!!!!!!!!!!\n");

        for ( int k = 0; k < data->n; k++ ) {
            if( data->tlk[k]->sm->m->freqs == stlk->sm->m->freqs ){
                data->tlk[k]->sm->m->set_frequency( data->tlk[k]->sm->m, Parameters_value(param, 0), i );
                SingleTreeLikelihood_update_all_nodes(data->tlk[k]);
            }
        }
        
        //lnl = fabs(stlk->calculate(stlk));
    }
	return -lnl;
}

// Make sure that one of the 2 branches below the root is set to 0 and is not optimized
static double _optimize_brent_branch_length_all( Optimizer *opt, BrentData2 *data, Parameters *param, int iterations ){
	opt_result status = OPT_NEED_CHECK;
    double lnl = 0;
    SingleTreeLikelihood *stlk = data->tlk[0];
    
    Node **nodes = Tree_get_nodes(stlk->tree, POSTORDER);
    for (int j = 0; j < iterations; j++){
        for (int i = 0; i < Tree_node_count(stlk->tree)-1; i++) {
            if( Parameter_fixed( nodes[i]->distance ) ){
                continue;
            }
            data->index_param = i;
            Parameters_set(param, 0, nodes[i]->distance);
            status = opt_optimize( opt, param, &lnl);
            if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
            Node_set_distance(nodes[i], Parameters_value(param, 0));
            
            for ( int k = 0; k < data->n; k++ ) {
                SingleTreeLikelihood_update_one_node( data->tlk[k], nodes[i]);
            }
        }
    }
	//return stlk->calculate(stlk);
    return -lnl;
}

double optimize_treelikelihoods( SingleTreeLikelihood **tlk, int nModels ){
	
    Tree *tree = tlk[0]->tree;
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	
	
	BrentData2 *data_brent = new_BrentData2( tlk, nModels );
	Parameters *oneparameter = new_Parameters(1);
    
    
	MultivariateData *data_multivariate = NULL;
    
	Optimizer *opt_freq = NULL;
	
	Optimizer *opt_rel_rate = NULL;
	
	Optimizer *opt_pinv = NULL;
	
	Optimizer *opt_gamma = NULL;
	
	Optimizer *opt_bl = NULL;
	
	OptConfig opt = tlk[0]->opt;
    
    Optimizer *opt_mu = new_Optimizer(OPT_BRENT);
    opt_set_data(opt_mu, data_brent );
    opt_set_objective_function(opt_mu, _brent_optimize_mu );
    
    
	
	// Init frequencies of the substitution model
	if( opt.freqs.optimize ){
		int n = Parameters_count_optimizable(tlk[0]->sm->m->freqs, NULL );
        for ( int i = 1; i < nModels; i++ ) {
            n += Parameters_count_optimizable(tlk[i]->sm->m->freqs, NULL );
        }
		
		if( n == 0 ){
			opt.freqs.optimize = false;
		}
		else{
			opt_freq = new_Optimizer( opt.freqs.method );
			opt_set_max_iteration(opt_freq, opt.freqs.max_iteration);
			
			if( opt.freqs.method == OPT_BRENT ){
				opt_set_data(opt_freq, data_brent );
				opt_set_objective_function(opt_freq, _brent_optimize_frequency );
				opt_set_tolx(opt_freq, opt.freqs.tolx);
			}
//			else if( opt.freqs.method == OPT_CG_PR || opt.freqs.method == OPT_CG_FR ){
//				data_multivariate = new_MultivariateData( stlk, NULL);
//				opt_set_data(opt_freq, data_multivariate);
//				opt_set_objective_function(opt_freq, _cg_optimize_frequencies);
//				opt_set_tolfx(opt_freq, opt.freqs.tolfx);
//			}
			else{
				fprintf(stderr, "Optmization not supported: frequencies %d %d\n", opt.freqs.method, OPT_BRENT  );
				exit(1);
			}
		}
	}
	
	// Init relative rates of the subsitution model
	if( opt.relative_rates.optimize ){
		int n = Parameters_count_optimizable(tlk[0]->sm->m->rates, NULL );
		for ( int i = 1; i < nModels; i++ ) {
            n += Parameters_count_optimizable(tlk[i]->sm->m->rates, NULL );
        }
        
		if( n == 0 ){
			opt.relative_rates.optimize = false;
		}
		else{
			opt_rel_rate = new_Optimizer( opt.relative_rates.method );
			opt_set_max_iteration(opt_rel_rate, opt.relative_rates.max_iteration );
			
			if( n == 1 || opt.relative_rates.method == OPT_BRENT ){
				opt_set_data(opt_rel_rate, data_brent );
				opt_set_objective_function(opt_rel_rate, _brent_optimize_relative_rate );
			}
            
//			else if( opt.relative_rates.method == OPT_CG_PR || opt.relative_rates.method == OPT_CG_FR ){
//				if( data_multivariate == NULL ) data_multivariate = new_MultivariateData( stlk, NULL);
//				opt_set_data(opt_rel_rate, data_multivariate);
//				opt_set_objective_function(opt_rel_rate, _cg_optimize_rel_rates);
//				opt_set_tolfx(opt_rel_rate, opt.relative_rates.tolfx);
//			}
			else{
				fprintf(stderr, "Optmization not supported: relative rates\n" );
				exit(1);
			}
		}
        
	}
	
	
	// Init proportion invariant
	if ( opt.pinv.optimize ) {
        int i = 0;
        for ( int i = 0; i < nModels; i++ ) {
            if( tlk[i]->sm->pinv != NULL ) break;
        }
        
		if( i != nModels ){
			opt_pinv = new_Optimizer( OPT_BRENT );
			opt_set_data(opt_pinv, data_brent );
			opt_set_objective_function(opt_pinv, _brent_optimize_pinv );
			//opt_set_max_iteration(opt_pinv, opt.pinv.max_iteration);
			//opt_set_tolx(opt_pinv, opt.pinv.tolx);
		}
		else {
			opt.pinv.optimize = false;
		}
        
	}
	
	// Init rate heterogeneity
	if ( opt.gamma.optimize ) {
        int i = 0;
        for ( int i = 0; i < nModels; i++ ) {
            if( tlk[i]->sm->shape != NULL ) break;
        }
        
		if( i != nModels ){
			opt_gamma = new_Optimizer( OPT_BRENT );
			opt_set_data(opt_gamma, data_brent);
			opt_set_objective_function(opt_gamma, _brent_optimize_gamma );
			//opt_set_max_iteration(opt_gamma, opt.gamma.max_iteration);
			//opt_set_tolx( opt_gamma, opt.gamma.tolx);
		}
		else {
			opt.gamma.optimize = false;
		}
	}
	
	
	// Init branch length optimization
	if( opt.bl.optimize ){
        // fix the right node, set it to 0 and add its distance to the left node
        // we should unfix it at the end
        Node *right = Node_right( Tree_root(tree) );
        Node *left  = Node_left( Tree_root(tree) );
        Node_set_distance(left, Node_distance(right)+Node_distance(left) );
        Node_set_distance(right, 0);
        Parameter_set_fixed(right->distance, true);
        
        for ( int i = 0; i < nModels; i++ ) {
            SingleTreeLikelihood_update_one_node(tlk[i], right);
            SingleTreeLikelihood_update_one_node(tlk[i], left);
        }
        
        int count = 0;
        for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
            if( !Parameter_fixed(nodes[i]->distance) ){
                count++;
            }
        }
		
		if( count == 0 ){
			opt.bl.optimize = false;
		}
		else{
            if(opt.bl.method == OPT_BRENT ){
                opt_bl= new_Optimizer( OPT_BRENT );
                opt_set_data(opt_bl, data_brent);
                opt_set_objective_function(opt_bl, _brent_optimize_branch_length);
                opt_set_max_iteration(opt_bl, opt.bl.max_iteration);
                opt_set_tolx( opt_bl, opt.bl.tolx);
            }
            else {
                error("Optmization not supported: branch length\n" );
            }
        }
        
        
	}
	
	int max_rounds = opt.max_rounds;
    
	//double lnl = stlk->calculate(stlk);
    double lnl = data_brent->f(data_brent);
	
	if( !opt.freqs.optimize && !opt.relative_rates.optimize && !opt.rates.optimize && !opt.heights.optimize && !opt.pinv.optimize && !opt.gamma.optimize && !opt.bl.optimize ){
		max_rounds = -1;
	}
	

	double fret = lnl;
	
	//opt.verbosity = 3;
    
	if ( opt.verbosity > -1 ) {
		fprintf(stderr, "Optimize branch lengths %s\n", (opt.bl.optimize ? "yes" : "no"));
		fprintf(stderr, "Optimize frequencies    %s\n", (opt.freqs.optimize ? "yes" : "no"));
		fprintf(stderr, "Optimize gamma          %s\n", (opt.gamma.optimize ? "yes" : "no"));
		fprintf(stderr, "Optimize pinv           %s\n", (opt.pinv.optimize ? "yes" : "no"));
		fprintf(stderr, "Optimize relative rates %s\n", (opt.relative_rates.optimize ? "yes" : "no"));
		
		fprintf(stderr, "Optimize heights        %s\n", (opt.heights.optimize ? "yes" : "no"));
		fprintf(stderr, "Optimize rates          %s\n", (opt.rates.optimize ? "yes" : "no"));
		fprintf(stderr, "Optimize topology       %s\n", (opt.topology_optimize ? "yes" : "no"));
		fprintf(stderr, "\n");
	}
	
	int iter_bl = 2;
	
    
    if ( opt.verbosity ==  1 ) {
        fprintf(stderr, "LnL: %f\n", lnl);
        
    }
    
	for ( int rounds = 0; rounds < max_rounds; rounds++ ) {
		//fret = 0;
		if ( opt.verbosity > 1  ) {
			fprintf(stderr, "========================================\n");
			fprintf(stderr, "Round %d lk = %f\n\n", rounds, lnl);
		}
		
		if ( opt.verbosity > 0 ) {
			fprintf(stderr, "\n");
		}
        
		// Branch length
		if( opt.bl.optimize ){
            function_f f = data_brent->f;
            if( tlk[0]->use_upper ){
                f = data_brent->f;
                data_brent->f = _standard_loglikelihood_upper_brent;
            }
            
			if( opt.bl.method == OPT_BRENT ){
				fret = _optimize_brent_branch_length_all(opt_bl, data_brent, oneparameter, iter_bl);
                
				if ( opt.verbosity > 0 ) {
					fprintf(stderr, "Branch lengths LnL: %f tree length %f {%f}\n", fret, Tree_length(tree), fret-lnl );
				}
				
				if ( fret - lnl < 0.1 ) {
					iter_bl = 1;
				}
				
			}
			else{
				error("Optimize: Only Brent is implemented for branch length optimization\n");
			}
            data_brent->f = f;
		}
        
		// Frequencies
		if( opt.freqs.optimize ){
            for ( int i = 0; i < nModels; i++) {
                
                int j = 0;
                for ( ; j < i; j++) {
                    if( tlk[i]->sm->m->freqs == tlk[j]->sm->m->freqs ){
                        break;
                    }
                }
                
                if ( j != i ) continue;
                
                SingleTreeLikelihood *stlk = tlk[i];
                
                if( opt.freqs.method == OPT_BRENT ){
                    fret = _optimize_brent_frequencies_all(i, opt_freq, data_brent, oneparameter);
                }
                else {
                    double status = opt_optimize( opt_freq, stlk->sm->m->freqs, &fret);
                    if( status == OPT_ERROR ) error("OPT.frequencies No SUCCESS!!!!!!!!!!!!\n");
                    fret = -fret;
                }
                
                if ( opt.verbosity > 0 ) {
                    stlk->sm->m->update_frequencies(stlk->sm->m);
                    fprintf(stderr, "Frequencies %d     LnL: %f [", i, fret);
                    for ( int i = 0; i < stlk->sm->nstate; i++ ) {
                        fprintf( stderr, "%s%f", (i==0?"":","), stlk->sm->m->_freqs[i] );
                    }
                    fprintf(stderr, "]  {%f}\n", fret-lnl );
                }
            }
		}
		
		// Relative rates
		if( opt.relative_rates.optimize ){
            for ( int i = 0; i < nModels; i++) {
                
                int j = 0;
                for ( ; j < i; j++) {
                    if( tlk[i]->sm->m->rates == tlk[j]->sm->m->rates ){
                        break;
                    }
                }
                
                if ( j != i ) continue;
                
                SingleTreeLikelihood *stlk = tlk[i];
                
                if( opt.relative_rates.method == OPT_BRENT ){
                    fret = _optimize_brent_relative_rate_all(i, opt_rel_rate, data_brent, oneparameter);
                }
                else{
                    double status = opt_optimize( opt_rel_rate, stlk->sm->m->rates, &fret);
                    if( status == OPT_ERROR ) error("OPT.relative_rates No SUCCESS!!!!!!!!!!!!\n");
                    fret = -fret;
                }
                
                if ( opt.verbosity > 0 ) {
                    fprintf(stderr, "Relative rates %d  LnL: %f [", i, fret);
                    for ( int i = 0; i < Parameters_count(stlk->sm->m->rates); i++ ) {
                        fprintf( stderr, "%s%f", (i==0?"":","), Parameters_value(stlk->sm->m->rates,i) );
                    }
                    fprintf(stderr, "]  {%f} [%d]\n", fret-lnl,opt.relative_rates.method );
                }
            }
		}
        
		// Gamma distributed rate heterogeneity
		if ( opt.gamma.optimize ) {
            
            for ( int i = 0; i < nModels; i++) {
                
                int j = 0;
                for ( ; j < i; j++) {
                    if( tlk[i]->sm->shape == tlk[j]->sm->shape ){
                        break;
                    }
                }
                if ( j != i ) continue;
                
                SingleTreeLikelihood *stlk = tlk[i];
                data_brent->index_model = i;
                
                // we ignore constraints here
                double lower = Parameter_lower(stlk->sm->shape);
                double upper = Parameter_upper(stlk->sm->shape);
                double value = Parameter_value(stlk->sm->shape);
                
                Parameter_set_lower(stlk->sm->shape, dmax(SITEMODEL_ALPHA_MIN, value/2) );
                Parameter_set_upper(stlk->sm->shape, dmin(SITEMODEL_ALPHA_MAX, value*2) );
                
                Parameters_set( oneparameter, 0,  stlk->sm->shape );
                
                double status= opt_optimize( opt_gamma, oneparameter, &fret);
                if( status == OPT_ERROR ) error("OPT.GAMMA No SUCCESS!!!!!!!!!!!!\n");
                
                for ( j = 0; j < nModels; j++ ) {
                    if( tlk[j]->sm->shape == stlk->sm->shape ){
                        SiteModel_set_alpha( tlk[j]->sm, Parameters_value(oneparameter, 0) );
                        SingleTreeLikelihood_update_all_nodes(tlk[j]);
                    }
                }
                
                fret = -fret;
                
                if ( opt.verbosity > 0 ) {
                    fprintf(stderr, "Gamma %d           LnL: %f shape: %f  {%f}\n", i, fret, Parameter_value(stlk->sm->shape), fret-lnl );
                }
                
                Parameter_set_bounds(stlk->sm->shape, lower, upper);
            }
		}
		
		// Proportion of invariant sites
		if ( opt.pinv.optimize ) {
            for ( int i = 0; i < nModels; i++) {
                
                int j = 0;
                for ( ; j < i; j++) {
                    if( tlk[i]->sm->pinv == tlk[j]->sm->pinv ){
                        break;
                    }
                }
                
                if ( j != i ) continue;
                
                SingleTreeLikelihood *stlk = tlk[i];
                data_brent->index_model = i;
                
                Parameters_set( oneparameter, 0,  stlk->sm->pinv );
                
                double status= opt_optimize( opt_pinv, oneparameter, &fret);
                if( status == OPT_ERROR ) error("OPT.PINV No SUCCESS!!!!!!!!!!!!\n");
                
                for ( j = 0; j < nModels; j++ ) {
                    if( tlk[j]->sm->pinv == stlk->sm->pinv ){
                        SiteModel_set_pinv( tlk[j]->sm, Parameters_value(oneparameter, 0) );
                        SingleTreeLikelihood_update_all_nodes(tlk[j]);
                    }
                }
                
                fret = -fret;
                
                if ( opt.verbosity > 0 ) {
                    fprintf(stderr, "Prop invariant %d  LnL: %f p: %f  {%f}\n", i, fret, Parameter_value(stlk->sm->pinv), fret-lnl );
                }
            }
		}
		
        if( nModels > 1 ){
            for ( int i = 0; i < nModels; i++) {
                SingleTreeLikelihood *stlk = tlk[i];
                data_brent->index_model = i;
                
                Parameters_set( oneparameter, 0,  stlk->sm->mu );
                
                double status= opt_optimize( opt_mu, oneparameter, &fret);
                if( status == OPT_ERROR ) error("OPT.MU No SUCCESS!!!!!!!!!!!!\n");
                
                
                SiteModel_set_mu( stlk->sm, Parameters_value(oneparameter, 0) );
                SingleTreeLikelihood_update_all_nodes(stlk);                  
                
                fret = -fret;
                
                if ( opt.verbosity > 0 ) {
                    fprintf(stderr, "Mu %d              LnL: %f [%f]  {%f}\n", i, fret, Parameter_value(stlk->sm->mu), fret-lnl );
                }
            }
        }
		
		if(  fret - lnl < opt.precision && rounds > 2 ){
            if( opt.verbosity > 1 ){
                fprintf(stderr, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
                fprintf(stderr, "STOP LnL previous %f last %f round: %d\n\n", lnl, fret, rounds);
            }
			lnl = fret;
			break;
		}
		
		lnl = fret;
        
	}
    
    free_BrentData2(data_brent);
    if (data_multivariate != NULL ) free_MultivariateData( data_multivariate);
    
    free_Parameters_soft(oneparameter);
    
	if( opt.freqs.optimize ){
		free_Optimizer( opt_freq );
	}
	
	if( opt.relative_rates.optimize ){
		free_Optimizer( opt_rel_rate );
	}
	
	if ( opt.pinv.optimize ) {
		free_Optimizer( opt_pinv );
	}
	
	if ( opt.gamma.optimize ) {
		free_Optimizer( opt_gamma );
	}
	
	if ( opt.bl.optimize ){
		free_Optimizer( opt_bl );
		
	}
    
    free_Optimizer(opt_mu);
    
	return lnl;
}
