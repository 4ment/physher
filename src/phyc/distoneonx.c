//
//  distoneonx.c
//  physher
//
//  Created by mathieu on 9/7/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distoneonx.h"

double DistributionModel_one_on_x_log_prob(DistributionModel* dm) {
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    Parameter* x = Parameters_at(dm->x, 0);
    for(size_t i = 0; i < Parameter_size(x); i++){
        dm->lp -= log(Parameter_value(x));   
    }
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_one_on_x_log_prob_grad(DistributionModel* dm, const Parameters* parameters) {
    Parameter* x = Parameters_at(dm->x, 0);
    Parameter* xx = Parameters_depends(parameters, x);
    if (xx != NULL) {
        double dlogP = -1.0 / Parameter_value(x);
        x->grad[0] += dlogP;
        if(xx != x){
            x->transform->backward(x->transform, &dlogP);
        }
    }
    return 0;
}

DistributionModel* new_OneOnXDistributionModel(Parameters* x) {
    DistributionModel* dm = new_DistributionModel(NULL, x);
    dm->type = DISTRIBUTION_ONE_ON_X;
    dm->log_prob = DistributionModel_one_on_x_log_prob;
    dm->log_prob_grad = DistributionModel_one_on_x_log_prob_grad;
    return dm;
}

Model* new_OneOnXDistributionModel_from_json(json_node* node, Hashtable* hash) {
    char* id = get_json_node_value_string(node, "id");

    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);

    DistributionModel* dm = new_OneOnXDistributionModel(x);

    Model* model = new_DistributionModel2(id, dm);

    model->samplable = false;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    free_Parameters(x);
    return model;
}
