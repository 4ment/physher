#include "minunit.h"
#include "phyc/filereader.h"
#include "phyc/hashtable.h"
#include "phyc/matrix.h"
#include "phyc/treelikelihood.h"

char* test_treelikelihood_time() {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);

    char* content = load_file("jc69-time.json");
    json_node* json = create_json_tree(content);
    free(content);

    json_node* child = json->children[0];
    Model* model = new_TreeLikelihoodModel_from_json(child, hash);
    SingleTreeLikelihood* tlk = model->obj;
    Model** models = (Model**)model->data;
    Model* mtree = models[0];
    Model* mbm = models[3];
    Tree* tree = mtree->obj;
    BranchModel* bm = mbm->obj;
    Node** nodes = Tree_nodes(tree);
    // init_heights_from_bls(tree);
    Tree_update_heights(tree);

    tlk->include_jacobian = false;

    double expected_logP = -4777.616349713985;
    double logP = model->logP(model);
    mu_assert(fabs(logP - expected_logP) < 1.e-8, "logP not matching");

    int flags = TREELIKELIHOOD_FLAG_TREE_MODEL | TREELIKELIHOOD_FLAG_BRANCH_MODEL;
    double* gradient = dvector(69);
    
    Parameters* params = get_reparams(tree);
    Parameter* ratios = Parameters_at(params, 0);
    Parameter* root = Parameters_at(params, 1);
    Parameters* parameters = new_Parameters(3);
    Parameters_add_parameters(parameters, params);
    Parameters_add_parameters(parameters, mbm->parameters);

    double expected_rate_grad = 328017.6732813406;
    double expected_ratio_grad[67] = {
        -0.5936536642214764, 6.441289658869611,   8.92145177998445,
        5.173924439035883,   -5.1189486033502325, 2.7314018967274634,
        2.007882472548766,   3.956031262797951,   5.542287760475186,
        9.56623809386586,    15.27690567000365,   35.18003581182256,
        73.00436877780763,   96.69564894572747,   14.99114774606325,
        15.285818508377771,  -1.3363345353505567, 10.94108984814406,
        19.64314696205841,   21.460133409615363,  39.1394523375063,
        3.637275922119337,   11.269174317983369,  12.443235860074363,
        71.12758013218424,   -3.8069961277876336, 88.1258829065779,
        3.5996001830340103,  18.479485706097613,  6.036534490720715,
        19.841103281559672,  23.24734623488343,   22.7331642319324,
        1.8172474126372273,  9.368306385819489,   54.08739297309535,
        42.35386071758409,   10.679777674119268,  4.140801615932186,
        3.3305556707250425,  -4.622247216603871,  27.32069418310099,
        54.31412932090593,   152.27137882559083,  23.540874887614432,
        14.3065705842615,    1.2225681560992132,  16.980030076368237,
        26.38017246149551,   3.4861149347888336,  4.098873332100652,
        10.267812216719863,  15.592298788222287,  70.94321518451146,
        4.240029132899654,   6.016353791291106,   38.343497684323275,
        3.4885156350078015,  66.51533636215693,   7.694985489230656,
        5.883423757661899,   3.981016102813299,   5.47007162703107,
        40.51912724901265,   30.451660702191045,  2.840830939900187,
        6.802521820384058};
    double expected_root_height_grad = 17.492484957839924;

// do it several times to check the gradient is properly reset
for(size_t k = 0 ; k < 2; k++){
    TreeLikelihood_gradient(model, flags, gradient);
    Parameters_zero_grad(parameters);
    model->gradient(model, parameters);

    mu_assert(fabs(gradient[68] - expected_rate_grad) < 1.e-8,
              "dlogP/dx clock not matching with TreeLikelihood_gradient");

    mu_assert(fabs(Parameters_at(parameters, 2)->grad[0] - expected_rate_grad) < 1.e-8,
              "dlogP/dx clock not matching");

    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(gradient[i] - expected_ratio_grad[i]) < 1.e-8,
                  "dlogP ratios not matching with TreeLikelihood_gradient");
        mu_assert(fabs(ratios->grad[i] - expected_ratio_grad[i]) < 1.e-8,
                  "dlogP ratios not matching");
    }

    mu_assert(fabs(gradient[67] - expected_root_height_grad) < 1.e-8,
              "dlogP root not matching with TreeLikelihood_gradient");
    
    mu_assert(fabs(root->grad[0] - expected_root_height_grad) < 1.e-8,
              "dlogP root not matching");
}

    // check the gradient is accumulated in parameter->grad
    // Parameters_zero_grad is not called so the gradient should be doubled
    model->gradient(model, parameters);

    mu_assert(fabs(Parameters_at(parameters, 2)->grad[0] - expected_rate_grad*2) < 1.e-8,
              "dlogP/dx clock not matching in accumulated");

    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(ratios->grad[i] - expected_ratio_grad[i]*2) < 1.e-8,
                  "dlogP ratios not matching in accumulated");
    }
    
    mu_assert(fabs(root->grad[0]  - expected_root_height_grad*2) < 1.e-8,
              "dlogP root not matching in accumulated");

    tlk->include_jacobian = true;

    double expected_logP_jacobian = -4786.867701371271;
    logP = model->logP(model);
    mu_assert(fabs(logP - expected_logP_jacobian) < 1.e-8, "logP not matching");

    double expected_ratios_jac_grad[67] = {
        -0.5936536642214764, 6.441289658869611,   11.202945298115116,
        5.173924439035883,   -0.9046311891428063, 2.7314018967274634,
        3.1571313705195485,  7.082913909386436,   10.305417331645046,
        13.988205820544293,  20.709336065224214,  48.897992914081215,
        99.16494936812502,   130.20574669099852,  17.314018642574176,
        21.033289555358838,  -1.3363345353505567, 12.259822362587805,
        22.88729131298567,   27.17656445923329,   47.48742627517851,
        3.637275922119337,   12.955169498485168,  15.31595344286499,
        83.25460505860441,   -3.8069961277876336, 105.38509458853852,
        4.874022850066035,   22.754466304821086,  6.036534490720715,
        25.651478211887106,  29.535185027483895,  29.598789450352278,
        1.8172474126372273,  10.598684711100873,  76.25924840292916,
        56.481422939218746,  10.679777674119268,  6.5871791334230085,
        3.3305556707250425,  -4.622247216603871,  33.41730442097831,
        63.4157671002785,    188.80951477041825,  23.540874887614432,
        17.42107593719064,   1.2225681560992132,  22.37201215315777,
        34.239511260483326,  3.4861149347888336,  4.098873332100652,
        13.200954262988732,  19.726890439483917,  96.80873776982577,
        4.240029132899654,   7.414584510049101,   48.87169351223057,
        3.4885156350078015,  82.96906517317022,   9.009333759946228,
        8.032474365132352,   3.981016102813299,   6.543650266561743,
        53.70242275637265,   37.835952010113665,  2.840830939900187,
        7.517186267961684};
    double expected_root_height_jac_grad = 19.936860572419484;
    SingleTreeLikelihood_update_all_nodes(tlk);
    double* gradientWithJacobian = dvector(69);
    TreeLikelihood_gradient(model, flags, gradientWithJacobian);

    mu_assert(fabs(gradientWithJacobian[68] - expected_rate_grad) < 1.e-8,
              "dlogP clock include_jacobian not matching");

    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(gradientWithJacobian[i] - expected_ratios_jac_grad[i]) < 1.e-8,
                  "dlogP ratios include_jacobian not matching");
    }

    mu_assert(fabs(gradientWithJacobian[67] - expected_root_height_jac_grad) < 1.e-8,
              "dlogP root include_jacobian not matching");

    // check log det jacobian and its gradient from tree Model alone
    double logDetJacobian = mtree->logP(mtree);
    double* logDetJacbianGrad = dvector(69);
    mu_assert(fabs(expected_logP_jacobian - expected_logP - logDetJacobian) < 1.e-8,
              "log det jacobian not matching");
    // Parameters* params = get_reparams(tree);
    // Parameter* ratios = Parameters_at(params, 0);
    // Parameter* root = Parameters_at(params, 1);
    Parameters_zero_grad(params);
    mtree->gradient(mtree, params);

    mu_assert(fabs(gradientWithJacobian[67] - gradient[67] - root->grad[0]) < 1.e-8,
              "dlogP root det jacobian not matching");
    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(gradientWithJacobian[i] - gradient[i] - ratios->grad[i]) < 1.e-8,
                  "dlogP ratios det jacobian not matching");
    }

    free_Parameters(parameters);
    model->free(model);
    // Model* sequenceModel = Hashtable_get(hash, "seqs");
    // sequenceModel->free(sequenceModel);
    free(gradient);
    free(gradientWithJacobian);
    free(logDetJacbianGrad);
    free_Hashtable(hash);
    json_free_tree(json);
    return NULL;
}

char* test_treelikelihood_time_unconstrained() {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);
    char* content = load_file("jc69-time.json");
    json_node* json = create_json_tree(content);
    free(content);
    json_node* child = json->children[1];
    Model* model = new_TreeLikelihoodModel_from_json(child, hash);
    SingleTreeLikelihood* tlk = model->obj;
    Model** models = (Model**)model->data;
    Model* mtree = models[0];
    Model* mbm = models[3];
    Tree* tree = mtree->obj;
    BranchModel* bm = mbm->obj;
    Node** nodes = Tree_nodes(tree);
    // init_heights_from_bls(tree);
    Tree_update_heights(tree);

    tlk->include_jacobian = false;

    double expected_logP = -4777.616209217359;
    double logP = model->logP(model);
    mu_assert(fabs(logP - expected_logP) < 1.e-8, "logP not matching");

    int flags = TREELIKELIHOOD_FLAG_TREE_MODEL | TREELIKELIHOOD_FLAG_BRANCH_MODEL;
    double* gradient = dvector(69);
    
    Parameters* params = get_reparams(tree);
    Parameter* ratios = Parameters_at(params, 0);
    Parameter* root = Parameters_at(params, 1);
    Parameter* rate = Parameters_at(mbm->parameters, 0);
    Parameters* parameters = new_Parameters(3);
    Parameters_add(parameters, ratios->transform->parameter);
    Parameters_add(parameters, root->transform->parameter);
    Parameters_add(parameters, rate->transform->parameter);
    // Parameters_add_parameters(parameters, mbm->parameters);

    double expected_rate_grad = 328017.6732813406;
    double expected_rate_unconstrained_grad = 328.01765057378356;
    double expected_ratio_grad[67] = {
        -0.5936536642214764, 6.441289658869611,   8.92145177998445,
        5.173924439035883,   -5.1189486033502325, 2.7314018967274634,
        2.007882472548766,   3.956031262797951,   5.542287760475186,
        9.56623809386586,    15.27690567000365,   35.18003581182256,
        73.00436877780763,   96.69564894572747,   14.99114774606325,
        15.285818508377771,  -1.3363345353505567, 10.94108984814406,
        19.64314696205841,   21.460133409615363,  39.1394523375063,
        3.637275922119337,   11.269174317983369,  12.443235860074363,
        71.12758013218424,   -3.8069961277876336, 88.1258829065779,
        3.5996001830340103,  18.479485706097613,  6.036534490720715,
        19.841103281559672,  23.24734623488343,   22.7331642319324,
        1.8172474126372273,  9.368306385819489,   54.08739297309535,
        42.35386071758409,   10.679777674119268,  4.140801615932186,
        3.3305556707250425,  -4.622247216603871,  27.32069418310099,
        54.31412932090593,   152.27137882559083,  23.540874887614432,
        14.3065705842615,    1.2225681560992132,  16.980030076368237,
        26.38017246149551,   3.4861149347888336,  4.098873332100652,
        10.267812216719863,  15.592298788222287,  70.94321518451146,
        4.240029132899654,   6.016353791291106,   38.343497684323275,
        3.4885156350078015,  66.51533636215693,   7.694985489230656,
        5.883423757661899,   3.981016102813299,   5.47007162703107,
        40.51912724901265,   30.451660702191045,  2.840830939900187,
        6.802521820384058};
    double expected_ratios_unconstrained_grad[67] = {
        -0.13274170070753702, 0.7757537767050319, 0.9649292279997707, 1.1327897228892874, -0.24708410677114834, 0.6151267488296555, 0.2268931184341151, 0.9118957036219864, 1.2921405461541102, 0.8257510376411261, 1.1192659367874378, 5.58923885747713, 18.126603385870844, 22.959722080018743, 2.6418528126199847, 3.8183876666762875, -0.1671482466858141, 2.0052767589800475, 4.644206371457391, 5.35183120408595, 2.4316470575593487, 0.24918727212011263, 2.7195717134513298, 2.631790338953837, 1.9087540051352347, -0.8602399958831338, 14.146078574968916, 0.43684000020042285, 4.61798903325333, 1.4165351345947172, 4.129274435446002, 3.575131223924172, 2.2670550015479596, 0.3480403066106321, 1.425694491228683, 13.085434563221968, 8.639569166459847, 1.8781871705264628, 1.0007380786793565, 0.6984737528022062, -0.6201082243176413, 6.022394850305062, 13.446244364474069, 37.9104400919227, 0.8361441423692472, 2.141575868027888, 0.17398749882036438, 2.251165586010986, 5.4825191893696035, 0.741393657507105, 1.0088653000435974, 2.227362577707621, 3.1046176816560265, 17.130539364674693, 0.7844883947046429, 1.2255051052496124, 5.630899806429098, 0.5106427284716466, 13.755725570710087, 1.9110158608426027, 1.1140807177193426, 0.7144352770901691, 0.34921331789395427, 9.208466923664183, 7.0775363562310565, 0.5197912298495507, 1.254968860544365
    };
    double expected_root_height_grad = 17.492484957839924;
    double expected_root_height_unconstrained_grad = 27.19186079725137;

// do it several times to check the gradient is properly reset
for(size_t k = 0 ; k < 2; k++){
    // TreeLikelihood_gradient(model, flags, gradient);
    Parameters_zero_grad(parameters);
    model->gradient(model, parameters);

    mu_assert(fabs(rate->transform->parameter->grad[0] - expected_rate_unconstrained_grad) < 1.e-8,
              "dlogP/dx clock not matching");

    for (size_t i = 0; i < 67; i++) {
        // mu_assert(fabs(gradient[i] - expected_ratio_grad[i]) < 1.e-8,
        //           "dlogP ratios not matching with TreeLikelihood_gradient");
        mu_assert(fabs(expected_ratios_unconstrained_grad[i] - ratios->transform->parameter->grad[i]) < 1.e-8,
                  "dlogP ratios not matching");
    }

    // mu_assert(fabs(gradient[67] - expected_root_height_grad) < 1.e-8,
    //           "dlogP root not matching with TreeLikelihood_gradient");
    
    mu_assert(fabs(expected_root_height_unconstrained_grad - root->transform->parameter->grad[0]) < 1.e-8,
              "dlogP root not matching");
}

    // check the gradient is accumulated in parameter->grad
    // Parameters_zero_grad is not called so the gradient should be doubled
    model->gradient(model, parameters);

    mu_assert(fabs(expected_rate_unconstrained_grad*2 - rate->transform->parameter->grad[0]) < 1.e-8,
              "dlogP/dx clock not matching in accumulated");

    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(expected_ratios_unconstrained_grad[i]*2 - ratios->transform->parameter->grad[i]) < 1.e-8,
                  "dlogP ratios not matching in accumulated");
    }
    
    mu_assert(fabs(expected_root_height_unconstrained_grad*2 - root->transform->parameter->grad[0]) < 1.e-8,
              "dlogP root not matching in accumulated");

    tlk->include_jacobian = true;

    double expected_logP_jacobian = -4786.867511318756;
    logP = model->logP(model);
    mu_assert(fabs(logP - expected_logP_jacobian) < 1.e-8, "logP not matching");

    double expected_ratios_jac_grad[67] = {
        -0.5936536642214764, 6.441289658869611,   11.202945298115116,
        5.173924439035883,   -0.9046311891428063, 2.7314018967274634,
        3.1571313705195485,  7.082913909386436,   10.305417331645046,
        13.988205820544293,  20.709336065224214,  48.897992914081215,
        99.16494936812502,   130.20574669099852,  17.314018642574176,
        21.033289555358838,  -1.3363345353505567, 12.259822362587805,
        22.88729131298567,   27.17656445923329,   47.48742627517851,
        3.637275922119337,   12.955169498485168,  15.31595344286499,
        83.25460505860441,   -3.8069961277876336, 105.38509458853852,
        4.874022850066035,   22.754466304821086,  6.036534490720715,
        25.651478211887106,  29.535185027483895,  29.598789450352278,
        1.8172474126372273,  10.598684711100873,  76.25924840292916,
        56.481422939218746,  10.679777674119268,  6.5871791334230085,
        3.3305556707250425,  -4.622247216603871,  33.41730442097831,
        63.4157671002785,    188.80951477041825,  23.540874887614432,
        17.42107593719064,   1.2225681560992132,  22.37201215315777,
        34.239511260483326,  3.4861149347888336,  4.098873332100652,
        13.200954262988732,  19.726890439483917,  96.80873776982577,
        4.240029132899654,   7.414584510049101,   48.87169351223057,
        3.4885156350078015,  82.96906517317022,   9.009333759946228,
        8.032474365132352,   3.981016102813299,   6.543650266561743,
        53.70242275637265,   37.835952010113665,  2.840830939900187,
        7.517186267961684};
    double expectedRatiosUnconstrainedJacobianGrad[67] = {
-0.13274170070753702, 0.7757537767050319, 1.21169172883093, 1.1327897228892874, -0.04366588628505873, 0.6151267488296555, 0.35675961670794615, 1.6326662598918036, 2.4026265317592084, 1.2074522227651623, 1.5172741708664017, 7.768683424270772, 24.62213868393231, 30.916466070863855, 3.0512066192138088, 5.254102963777162, -0.1671482466858141, 2.2469733002806693, 5.411215635665562, 6.777422270224085, 2.9502877952728577, 0.24918727212011263, 3.1264502312284317, 3.2393806581383786, 2.2341904671223562, -0.8602399958831338, 16.91654919457692, 0.5915046719807964, 5.686298796323111, 1.4165351345947172, 5.338512731751605, 4.542116335273242, 2.9517261384675835, 0.3480403066106321, 1.6129368145500587, 18.44949801237873, 11.521383627363338, 1.8781871705264628, 1.5919706446687631, 0.6984737528022062, -0.6201082243176413, 7.366292245205609, 15.699485878822047, 47.00720431253316, 0.8361441423692472, 2.6077916454436796, 0.17398749882036438, 2.966019246405346, 7.115904052991468, 0.741393657507105, 1.0088653000435974, 2.8636392779338156, 3.927865433352233, 23.37624470106702, 0.7844883947046429, 1.5103167096490036, 7.177008637398779, 0.5106427284716466, 17.158445964147663, 2.237431941450079, 1.5210238338970705, 0.7144352770901691, 0.41775132843043744, 12.204533158328786, 8.793784454079837, 0.5197912298495507, 1.3868141466982251
    };
    double expectedRateUnconstrainedJacobianGrad = 328.01765057378356;
    double expected_root_height_jac_grad = 19.936860572419484;
    double expectedRootHeightUnconstrainedJacobianGrad = 30.99161526317372;

    SingleTreeLikelihood_update_all_nodes(tlk);

    Parameters_zero_grad(parameters);
    model->gradient(model, parameters);

    mu_assert(fabs(expectedRateUnconstrainedJacobianGrad - rate->transform->parameter->grad[0]) < 1.e-8,
              "dlogP clock include_jacobian not matching");

    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(ratios->transform->parameter->grad[i] - expectedRatiosUnconstrainedJacobianGrad[i]) < 1.e-8,
                  "dlogP ratios include_jacobian not matching");
    }

    mu_assert(fabs(root->transform->parameter->grad[0] - expectedRootHeightUnconstrainedJacobianGrad) < 1.e-8,
              "dlogP root include_jacobian not matching");

    // check log det jacobian and its gradient from tree Model alone
    double logDetJacobian = mtree->logP(mtree);
    double* logDetJacbianGrad = dvector(69);
    mu_assert(fabs(expected_logP_jacobian - expected_logP - logDetJacobian) < 1.e-8,
              "log det jacobian not matching");
    // Parameters* params = get_reparams(tree);
    // Parameter* ratios = Parameters_at(params, 0);
    // Parameter* root = Parameters_at(params, 1);
    Parameters_zero_grad(parameters);
    mtree->gradient(mtree, parameters);

    mu_assert(fabs(expectedRootHeightUnconstrainedJacobianGrad - expected_root_height_unconstrained_grad - root->transform->parameter->grad[0]) < 1.e-8,
              "dlogP root det jacobian not matching");
    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(expectedRatiosUnconstrainedJacobianGrad[i] - expected_ratios_unconstrained_grad[i] - ratios->transform->parameter->grad[i]) < 1.e-8,
                  "dlogP ratios det jacobian not matching");
    }

    // rate is independent from the tree model transformation so the gradient should be zero
    mu_assert(rate->transform->parameter->grad[0] == 0.0, "dlogP rate det jacobian not matching");

    free_Parameters(parameters);
    model->free(model);
    // Model* sequenceModel = Hashtable_get(hash, "seqs");
    // sequenceModel->free(sequenceModel);
    free(gradient);
    free(logDetJacbianGrad);
    free_Hashtable(hash);
    json_free_tree(json);
    return NULL;
}

// #include "phyc/distnormal.h"
// char* test_normal_distribution_issigma(bool issigma) {
    
//     Parameter* mu = new_Parameter("mu", 2.0, new_Constraint(-INFINITY, INFINITY));
//     Parameter* sigma = new_Parameter("sigma", 0.1, new_Constraint(0, INFINITY));
//     Parameters* parameters = new_Parameters(2);
//     Parameter* x = new_Parameter("x", 0.1, new_Constraint(0, INFINITY));
//     Parameters_add(ps, mu);
//     Parameters_add(ps, sigma);
//     DistributionModel* dm = NULL;
//     if (issigma) {
//         dm = new_NormalDistributionModel_with_parameters(
//             parameters, xs, DISTRIBUTION_NORMAL_MEAN_SIGMA);
//     } else {
//         dm = new_NormalDistributionModel_with_parameters(ps, xs,
//                                                          DISTRIBUTION_NORMAL_MEAN_TAU);
//     }
//     Model* model = new_DistributionModel2("dist", dm);

//     double logP = model->logP(model);
//     double logP2 = gsl_normal_logP(x, xdim, p, pdim, issigma);
//     mu_assert(logP == logP2, "logP not matching");

//     x[0] = 10;
//     Parameters_set_value(xs, 0, x[0]);
//     logP = model->logP(model);
//     logP2 = gsl_normal_logP(x, xdim, p, pdim, issigma);
//     mu_assert(logP == logP2, "logP after reset not matching");

//     double eps = 0.00001;
//     size_t idx = 1;
//     double dlogPdx = model->dlogP(model, Parameters_at(xs, idx));
//     double dlogPdx2 = gsl_normal_dlogPdx(x, xdim, idx, p, pdim, eps, issigma);
//     mu_assert(fabs(dlogPdx - dlogPdx2) < 0.0001, "dlogPdx not matching");

//     double d2logPdx = model->d2logP(model, Parameters_at(xs, idx));
//     double d2logPdx2 = gsl_normal_d2logPdx(x, xdim, idx, p, pdim, eps, issigma);
//     mu_assert(fabs(d2logPdx - d2logPdx2) < 0.0001, "d2logPdx not matching");
//     //
//     // 	double dlogPdmu = model->dlogP(model, Parameters_at(ps[0], 0));
//     // 	double dlogPdmu2 = gsl_normal_1dist_dlogPdp(xs, p, 0, eps, issigma);
//     // 	mu_assert(fabs(dlogPdmu - dlogPdmu2) < 0.0001, "dlogPdmu not matching");
//     //
//     // 	double dlogPds = model->dlogP(model, Parameters_at(ps[1], 0));
//     // 	double dlogPds2 = gsl_normal_1dist_dlogPdp(xs, p, 1, eps, issigma);
//     // 	mu_assert(fabs(dlogPds - dlogPds2) < 0.0001, "dlogPdsigma not matching");

//     model->free(model);
//     for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
//     free(ps);
//     return NULL;
// }

char* all_tests() {
    mu_suite_start();
    // mu_run_test(test_treelikelihood_time);
    mu_run_test(test_treelikelihood_time_unconstrained);

    return NULL;
}

RUN_TESTS(all_tests);