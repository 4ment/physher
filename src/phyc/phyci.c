/*
 *  phyci.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 14/11/2013.
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

#include "treelikelihood.h"
#include "matrix.h"
#include "chisq.h"
#include "optimize.h"
#include "optimizer.h"


static double _brent_optimize_rate_ci( Parameters *params, double *grad, void *data ){
	BrentData *mydata = (BrentData*)data;
	int index = mydata->index_param;
	SingleTreeLikelihood *stlk = mydata->tlk;
	
	stlk->bm->set( stlk->bm, index, Parameters_value(params, 0) );
    
    SingleTreeLikelihood_update_all_nodes(stlk);

    return fabs(mydata->f(mydata)-mydata->backup[0]);
}

double ** ConfidenceInterval_rates( SingleTreeLikelihood *tlk, double level ){
    double **ci = dmatrix(Parameters_count(tlk->bm->rates), 2 );
    double c = pchisq(1.0-level, 1);
    
    BrentData *data = new_BrentData(tlk);
    data->backup = dvector(1);
    data->backup[0] = tlk->calculate(tlk)-0.5*c;
    Optimizer *opt = new_Optimizer(OPT_BRENT);
    opt_set_data(opt, data);
    opt_set_objective_function(opt, _brent_optimize_rate_ci);
    opt_set_tolx(opt, 0.000000001);
    Parameters *param = new_Parameters(1);
    Parameters_add(param, new_Parameter("ci", 0, new_Constraint(0, 1)));
    double fmin = 0;
    for ( int i = 0; i < Parameters_count(tlk->bm->rates); i++ ) {
        data->index_param = i;
        double rate = Parameters_value(tlk->bm->rates, i);
        
        // lower
        Parameters_set_bounds(param, 0, rate/10, rate );
        Parameters_set_value(param, 0, rate);
        opt_optimize(opt, param, &fmin);
        ci[i][0] = Parameters_value(param, 0);
        
        Parameters_set_value(tlk->bm->rates, i, rate);
        
        // upper
        Parameters_set_bounds(param, 0, rate, rate*10);
        Parameters_set_value(param, 0, rate);
        opt_optimize(opt, param, &fmin);
        ci[i][1] = Parameters_value(param, 0);
        
        Parameters_set_value(tlk->bm->rates, i, rate);
    }
    free_Parameters(param);
    free_Optimizer(opt);
    free_BrentData(data);
    
    SingleTreeLikelihood_update_all_nodes(tlk);
    tlk->calculate(tlk);
    
    return ci;
}

