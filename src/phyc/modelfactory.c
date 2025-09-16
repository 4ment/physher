// #include "modelfactory.h"

// #include <string.h>

// #include "sequence.h"
// #include "demographicmodels.h"
// #include "compoundmodel.h"
// #include "distmodelfactory.h"
// #include "elbo.h"
// #include "klpq.h"
// #include "jacobiantransform.h"
// #include "parsimony.h"
// #include "tree.h"
// #include "treelikelihood.h"
// #include "mjson.h"

// Model* model_factory(json_node* node, Hashtable* hash){
//     const char* type = get_json_node_string(node, "type");
//     Model* model = NULL;

//     if (strcasecmp(type, "alignment") == 0) {
//         model = new_Alignment_from_json(node, hash);
//     }
//     else if (strcasecmp(type, "coalescent") == 0) {
//         model = new_CoalescentModel_from_json(node, hash);
//     }
//     else if (strcasecmp(type, "compound") == 0) {
//         model = new_CompoundModel_from_json(node, hash);
//     }
//     else if (strcasecmp(type, "distribution") == 0) {
//         model = new_DistributionModel_from_json(node, hash);
//     }
//     // else if(model_type == MODEL_VARIATIONAL){
//     // 	models[index] = new_Variational_from_json(child, hash2);
//     // }
//     else if (strcasecmp(type, "elbo") == 0) {
//         model = new_ELBO_from_json(node, hash);
//     }
//     else if (strcasecmp(type, "jacobian") == 0) {
//         model = new_JacobianTransformModel_from_json(node, hash);
//     }
//     else if (strcasecmp(type, "klpq") == 0) {
//         model = new_KLpq_from_json(node, hash);
//     }
//     else if (strcasecmp(type, "parsimony") == 0) {
//         model = new_ParsimonyModel_from_json(node, hash);
//     }
//     else if (strcasecmp(type, "tree") == 0) {
//         model = new_TreeModel_from_json(node, hash);
//     }
//     else if (strcasecmp(type, "treelikelihood") == 0) {
//         model = new_TreeLikelihoodModel_from_json(node, hash);
//     }
//     return model;
// }
