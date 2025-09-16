#include "distmodelfactory.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "ctmcscale.h"
#include "distbeta.h"
#include "distbetaprime.h"
#include "distcauchy.h"
#include "distdirichlet.h"
#include "distexp.h"
#include "distgamma.h"
#include "distkumaraswamy.h"
#include "distlognormal.h"
#include "distmultinormal.h"
#include "distnormal.h"
#include "distoneonx.h"
#include "exponential.h"
#include "filereader.h"
#include "gamma.h"
#include "gmrf.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"
#include "tree.h"
#include "utilsio.h"

Model* new_DistributionModel_from_json(json_node* node, Hashtable* hash) {
    char* allowed[] = {"burnin",     "distribution",     "file",      "from", "margin",
                       "parameters", "parameterization", "posterior", "tree", "x"};
    json_check_allowed(node, allowed, sizeof(allowed) / sizeof(allowed[0]));

    char* d_string = get_json_node_value_string(node, "distribution");

    if (strcasecmp(d_string, "beta") == 0) {
        return new_BetaDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "betaprime") == 0) {
        return new_BetaPrimeDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "exponential") == 0) {
        return new_ExponentialDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "gamma") == 0) {
        return new_GammaDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "gmrf") == 0) {
        return new_GMRFModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "dirichlet") == 0) {
        return new_DirichletDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "cauchy") == 0) {
        return new_CauchyDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "kumaraswamy") == 0) {
        return new_KumaraswamyDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "lognormal") == 0) {
        return new_LogNormalDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "multivariatenormal") == 0) {
        return new_MultivariateNormalDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "normal") == 0 ||
               strcasecmp(d_string, "gaussian") == 0) {
        return new_NormalDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "halfnormal") == 0) {
        return new_HalfNormalDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "oneonx") == 0) {
        return new_OneOnXDistributionModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "ctmcscale") == 0) {
        return new_CTMCScaleModel_from_json(node, hash);
    } else if (strcasecmp(d_string, "copula") == 0) {
        // return new_GaussianCopuaModel_from_json(node, hash);
    }

    char* id = get_json_node_value_string(node, "id");

    Parameters** parameters = NULL;
    size_t parameters_dim = 0;
    Parameters* x = new_Parameters(1);
    DistributionModel* dm = NULL;
    Model* model = NULL;

    char* file = get_json_node_value_string(node, "file");
    Vector** samples = NULL;
    if (file != NULL) {
        get_parameters_references2(node, hash, x, "x");
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        samples = read_log_for_parameters_t(file, burnin, x);
    }

    if (strcasecmp(d_string, "topology") == 0) {
        char* ref = get_json_node_value_string(node, "tree");
        Model* mtree = Hashtable_get(hash, ref + 1);
        dm = new_UniformTreeDistribution(mtree->obj);
        model = new_DistributionModel3(id, dm, mtree);
    } else if (strcasecmp(d_string, "oneonx") == 0) {
        get_parameters_references2(node, hash, x, "x");
        dm = new_OneOnXDistributionModel(x);
        model = new_DistributionModel2(id, dm);
    } else {
        printf("Distribution unknown: %s\n", d_string);
        exit(10);
    }

    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");

    if (samples != NULL) {
        size_t paramCount = Parameters_count(x);
        for (int i = 0; i < paramCount; i++) {
            free_Vector(samples[i]);
        }
        free(samples);
    }
    if (parameters != NULL) {
        for (int i = 0; i < parameters_dim; i++) free_Parameters(parameters[i]);
        free(parameters);
    }
    free_Parameters(x);

    return model;
}