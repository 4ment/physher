//
//  distoneonx.c
//  physher
//
//  Created by mathieu on 9/7/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distoneonx.h"

double DistributionModel_log_one_on_x(DistributionModel* dm) {
    return -log(Parameters_value(dm->x, 0));
}

double DistributionModel_log_one_on_x_with_values(DistributionModel* dm, const double* values) {
    return -log(values[0]);
}

double DistributionModel_dlog_one_on_x(DistributionModel* dm, const Parameter* p) {
    if (p == Parameters_at(dm->x, 0)) {
        return -1.0 / Parameters_value(dm->x, 0);
    }
    return 0;
}

double DistributionModel_d2log_one_on_x(DistributionModel* dm, const Parameter* p) {
    if (strcmp(Parameter_name(p), Parameters_name(dm->x, 0)) == 0) {
        return 1.0 / Parameters_value(dm->x, 0) / Parameters_value(dm->x, 0);
    }
    return 0;
}

DistributionModel* new_OneOnXDistributionModel(Parameters* x) {
    DistributionModel* dm = new_DistributionModel(NULL, 0, x);
    dm->type = DISTRIBUTION_ONE_ON_X;
    dm->logP = DistributionModel_log_one_on_x;
    dm->logP_with_values = DistributionModel_log_one_on_x_with_values;
    dm->dlogP = DistributionModel_dlog_one_on_x;
    dm->d2logP = DistributionModel_d2log_one_on_x;
    dm->ddlogP = DistributionModel_ddlog_0;
    return dm;
}

Model* new_OneOnXDistributionModel_from_json(json_node* node, Hashtable* hash) {
    char* id = get_json_node_value_string(node, "id");

    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);

    DistributionModel* dm = new_OneOnXDistributionModel(x);

    Model* model = new_DistributionModel2(id, dm);

    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    return model;
}
