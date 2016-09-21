/*
 *  mixedtreelikelihood.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/12/2015.
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



#include "mixedtreelikelihood.h"

#include <stdlib.h>
#include <assert.h>

#include "gamma.h"
#include "matrix.h"

#include "optimize.h"

static void _gamma_approx_quantile( Parameter *shape, int cat, double *rates, double *proportions );

static double _get_category_gamma( MixedTreeLikelihood *mixed, const int index );

static double _get_proportion_gamma( MixedTreeLikelihood *mixed, const int index );

static double _calculate( MixedTreeLikelihood *mixed );


MixedTreeLikelihood * new_MixedStrictTreeLikelihood( SingleTreeLikelihood *tlk, int n ){
    MixedTreeLikelihood *mixed = malloc(sizeof(MixedTreeLikelihood));
    assert(mixed);
    mixed->tlks = (SingleTreeLikelihood**)malloc(n * sizeof(SingleTreeLikelihood*) );
    assert(mixed->tlks);
    mixed->n = n;
    mixed->tlks[0] = tlk;
    
    for ( int i = 1; i < n; i++ ) {
        mixed->tlks[i] = new_SingleTreeLikelihood(tlk->tree, tlk->sm, tlk->sp, clone_BranchModel(tlk->bm, tlk->tree));
    }
    
    
    mixed->categories  = dvector(n);
    mixed->proportions = dvector(n);
    mixed->get_category   = _get_category_gamma;
    mixed->get_proportion = _get_proportion_gamma;
    mixed->lnl = 0;
    mixed->calculate = _calculate;
    mixed->params = new_Parameters(2);
    Parameter *alpha = new_Parameter("mixed.alpha", 1, new_Constraint(0.01, 100) );
	Parameters_add(mixed->params, alpha);
    
    double r = tlk->bm->get(tlk->bm, NULL);
    Parameter *rate = new_Parameter("mixed.rate", r, new_Constraint(r*0.01, r*20) );
	Parameters_add(mixed->params, rate);
    
    return mixed;
}

void free_MixedTreeLikelihood( MixedTreeLikelihood *mixed ){
    free(mixed->categories);
    free(mixed->proportions);
    for ( int i = 1; i < mixed->n; i++ ) {
        free_SingleTreeLikelihood_share2( mixed->tlks[i], true, true, true, false );
    }
    free_Parameters(mixed->params);
    free(mixed);
}

double _get_category_gamma( MixedTreeLikelihood *mixed, const int index ){
	if ( mixed->need_update ) {
		_gamma_approx_quantile(Parameters_at(mixed->params, 0), mixed->n, mixed->categories, mixed->proportions);
        mixed->need_update = false;
	}
	return mixed->categories[index];
}

double _get_proportion_gamma( MixedTreeLikelihood *mixed, const int index ){
	if ( mixed->need_update ) {
		_gamma_approx_quantile(Parameters_at(mixed->params, 0), mixed->n, mixed->categories, mixed->proportions);
        mixed->need_update = false;
	}
	return mixed->proportions[index];
}


void _gamma_approx_quantile( Parameter *shape, int cat, double *rates, double *proportions ) {
	double mean = 0.0;
    double p = 1.0 / cat;
    
    for (int i = 0; i < cat; i++) {
        rates[i] = qgamma( (2.0 * i + 1.0) / (2.0 * cat), Parameter_value(shape), 1.0 / Parameter_value(shape) );
        mean += rates[i];
        
        proportions[i] = p;
    }
    
    mean = mean/ cat;
    
    for (int i = 0; i < cat; i++) {
        rates[i] /= mean;
    }
}

double _calculate( MixedTreeLikelihood *mixed ){
    if( mixed->need_update ){
        _gamma_approx_quantile(Parameters_at(mixed->params, 0), mixed->n, mixed->categories, mixed->proportions);
        
        for ( int i = 0; i < mixed->n; i++ ) {
            mixed->tlks[i]->bm->set( mixed->tlks[i]->bm, 0, mixed->categories[i]*Parameters_value(mixed->params, 1) );
            SingleTreeLikelihood_update_all_nodes(mixed->tlks[i]);
        }
        mixed->need_update = false;
    }
    
    mixed->lnl = 0;
    for ( int i = 0; i < mixed->n; i++ ) {
        mixed->lnl += mixed->proportions[i] * mixed->tlks[i]->calculate(mixed->tlks[i]);
    }
    return mixed->lnl;
}


typedef struct GenericBrentData{
	void *tlk;
    function_f f;
	int index_param;
	double *backup;
} GenericBrentData;

double mixed_loglikelihood_brent( void *data ){
    MixedTreeLikelihood *mixed = (MixedTreeLikelihood*)((GenericBrentData*)data)->tlk;
    return mixed->calculate(mixed);
}

double _brent_optimize_params( Parameters *params, double *grad, void *data ){
	GenericBrentData *mydata = (GenericBrentData*)data;
	MixedTreeLikelihood *mixed = (MixedTreeLikelihood*)mydata->tlk;
	
    Parameters_set_value(mixed->params, mydata->index_param, Parameters_value(params, 0) );
    
    mixed->need_update = true;
    
    return fabs(mydata->f(mydata));
}

GenericBrentData * new_GenericBrentData( void *tlk ){
	GenericBrentData *data = (GenericBrentData*)malloc( sizeof(GenericBrentData) );
    assert(data);
	data->index_param = 0;
	data->tlk = tlk;
	data->backup = NULL;
    data->f = mixed_loglikelihood_brent;
    
	return data;
}

void free_GenericBrentData( GenericBrentData *data ){
	if( data->backup != NULL ){
		free( data->backup);
	}
	free(data);
}

double optimize_mixedtreelikelihoods( MixedTreeLikelihood *mixed ){
    GenericBrentData *data_brent = new_GenericBrentData( mixed );
    Optimizer *opt_brent = new_Optimizer( OPT_BRENT );
    opt_set_data(opt_brent, data_brent );
    opt_set_objective_function(opt_brent, _brent_optimize_params );
    opt_set_tolx(opt_brent, 0.001);
    
	Parameters *oneparameter = new_Parameters(1);
    
    double lnl = data_brent->f(data_brent);
    double fret;
    
    mixed->tlks[0]->opt.rates.optimize = false;
    
    for ( int rounds = 0; rounds < 500; rounds++ ) {
        
        // Gamma shape
        Parameters_set( oneparameter, 0,  Parameters_at(mixed->params, 0) );
        data_brent->index_param = 0;
        double status = opt_optimize( opt_brent, oneparameter, &fret);
        if( status == OPT_ERROR ) error("OPT No SUCCESS!!!!!!!!!!!!\n");
        fret = -fret;
        
        if ( mixed->tlks[0]->opt.verbosity > 0 ) {
            fprintf(stdout, "Gamma         LnL: %f shape: %f  {%f}\n", fret, Parameters_value(mixed->params,0), fret-lnl );
        }
        
        // Rate
        Parameters_set( oneparameter, 0,  Parameters_at(mixed->params, 1) );
        data_brent->index_param = 1;
        status = opt_optimize( opt_brent, oneparameter, &fret);
        if( status == OPT_ERROR ) error("OPT.rate No SUCCESS!!!!!!!!!!!!\n");
        fret = -fret;
        
        if ( mixed->tlks[0]->opt.verbosity > 0 ) {
            fprintf(stdout, "Rate          LnL: %f rate: %f  {%f}\n", fret, Parameters_value(mixed->params,1), fret-lnl );
        }
        
//        // Branches
//        fret = optimize_singletreelikelihood(mixed->tlks[0]);
//        if ( mixed->tlks[0]->opt.verbosity > 0 ) {
//            if ( mixed->tlks[0]->opt.verbosity > 1 ) fprintf(stderr, "---------------------------------------------------------------------\n");
//            fprintf(stdout, "Root height   LnL: %f height: %f {%f}\n", fret, Node_height(Tree_root(mixed->tlks[0]->tree)), fret-lnl );
//        }
        
        if(  fret - lnl < mixed->tlks[0]->opt.precision ){
            lnl = fret;
            break;
        }
        
    }
    free_Parameters_soft(oneparameter);
    free_GenericBrentData(data_brent);
    free_Optimizer( opt_brent );
    
    return lnl;
}

