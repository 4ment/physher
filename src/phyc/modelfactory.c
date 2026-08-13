// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "modelfactory.h"

#include <string.h>
#include <strings.h>

#include "model.h"
#include "sequence.h"
#include "boundfactory.h"
#include "demographicmodels.h"
#include "compoundmodel.h"
#include "distmodelfactory.h"
#include "jacobian.h"
#include "parsimony.h"
#include "tree.h"
#include "treelikelihood.h"
#include "mjson.h"
#include "sitemodel.h"
#include "substmodel.h"

Model* model_factory_from_json(json_node* node, Hashtable* hash){
    const char* type = get_json_node_value_string(node, "type");
    Model* model = NULL;

    if (strcasecmp(type, "alignment") == 0) {
        model = new_Alignment_from_json(node, hash);
    }
    else if (strcasecmp(type, "bound") == 0) {
        model = new_BoundModel_from_json(node, hash);
    }
    else if (strcasecmp(type, "coalescent") == 0) {
        model = new_CoalescentModel_from_json(node, hash);
    }
    else if (strcasecmp(type, "compound") == 0) {
        model = new_CompoundModel_from_json(node, hash);
    }
    else if (strcasecmp(type, "distribution") == 0) {
        model = new_DistributionModel_from_json(node, hash);
    }
    else if (strcasecmp(type, "jacobian") == 0) {
        model = new_JacobianModel_from_json(node, hash);
    }
    else if (strcasecmp(type, "parsimony") == 0) {
        model = new_ParsimonyModel_from_json(node, hash);
    }
    else if (strcasecmp(type, "sitemodel") == 0) {
        model = new_SiteModel_from_json(node, hash);
    }
    else if (strcasecmp(type, "substitutionmodel") == 0) {
        model = new_SubstitutionModel_from_json(node, hash);
    }
    else if (strcasecmp(type, "tree") == 0) {
        model = new_TreeModel_from_json(node, hash);
    }
    else if (strcasecmp(type, "treelikelihood") == 0) {
        model = new_TreeLikelihoodModel_from_json(node, hash);
    }
    return model;
}
