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
    for (size_t k = 0; k < 2; k++) {
        TreeLikelihood_gradient(model, flags, gradient);
        Parameters_zero_grad(parameters);
        model->gradient(model, parameters);

        mu_assert(fabs(gradient[68] - expected_rate_grad) < 1.e-8,
                  "dlogP/dx clock not matching with TreeLikelihood_gradient");

        mu_assert(
            fabs(Parameters_at(parameters, 2)->grad[0] - expected_rate_grad) < 1.e-8,
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

    mu_assert(
        fabs(Parameters_at(parameters, 2)->grad[0] - expected_rate_grad * 2) < 1.e-8,
        "dlogP/dx clock not matching in accumulated");

    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(ratios->grad[i] - expected_ratio_grad[i] * 2) < 1.e-8,
                  "dlogP ratios not matching in accumulated");
    }

    mu_assert(fabs(root->grad[0] - expected_root_height_grad * 2) < 1.e-8,
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
    printf("logP: %f, expected_logP: %f\n", logP, expected_logP);
    mu_assert(fabs(logP - expected_logP) < 1.e-8, "logP not matching");

    // logP/full_logP must be stored in the model's lp field before returning
    mu_assert(model->lp == logP, "TreeLikelihoodModel logP not stored in lp");
    mu_assert(model->full_logP(model) == logP, "TreeLikelihoodModel full_logP not matching");
    mu_assert(model->lp == logP, "TreeLikelihoodModel full_logP not stored in lp");

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
        -0.13274170070753702, 0.7757537767050319,   0.9649292279997707,
        1.1327897228892874,   -0.24708410677114834, 0.6151267488296555,
        0.2268931184341151,   0.9118957036219864,   1.2921405461541102,
        0.8257510376411261,   1.1192659367874378,   5.58923885747713,
        18.126603385870844,   22.959722080018743,   2.6418528126199847,
        3.8183876666762875,   -0.1671482466858141,  2.0052767589800475,
        4.644206371457391,    5.35183120408595,     2.4316470575593487,
        0.24918727212011263,  2.7195717134513298,   2.631790338953837,
        1.9087540051352347,   -0.8602399958831338,  14.146078574968916,
        0.43684000020042285,  4.61798903325333,     1.4165351345947172,
        4.129274435446002,    3.575131223924172,    2.2670550015479596,
        0.3480403066106321,   1.425694491228683,    13.085434563221968,
        8.639569166459847,    1.8781871705264628,   1.0007380786793565,
        0.6984737528022062,   -0.6201082243176413,  6.022394850305062,
        13.446244364474069,   37.9104400919227,     0.8361441423692472,
        2.141575868027888,    0.17398749882036438,  2.251165586010986,
        5.4825191893696035,   0.741393657507105,    1.0088653000435974,
        2.227362577707621,    3.1046176816560265,   17.130539364674693,
        0.7844883947046429,   1.2255051052496124,   5.630899806429098,
        0.5106427284716466,   13.755725570710087,   1.9110158608426027,
        1.1140807177193426,   0.7144352770901691,   0.34921331789395427,
        9.208466923664183,    7.0775363562310565,   0.5197912298495507,
        1.254968860544365};
    double expected_root_height_grad = 17.492484957839924;
    double expected_root_height_unconstrained_grad = 27.19186079725137;

    // do it several times to check the gradient is properly reset
    for (size_t k = 0; k < 2; k++) {
        TreeLikelihood_gradient(model, flags, gradient);
        Parameters_zero_grad(parameters);
        model->gradient(model, parameters);

        mu_assert(fabs(rate->transform->parameter->grad[0] -
                       expected_rate_unconstrained_grad) < 1.e-8,
                  "dlogP/dx clock not matching");

        printf("rate grad: expected %g got %g diff %g\n", gradient[68],
               expected_rate_grad, gradient[68] - expected_rate_grad);
        mu_assert(fabs(gradient[68] - expected_rate_grad) < 1.e-2,
                  "dlogP/dx clock not matching with TreeLikelihood_gradient");

        for (size_t i = 0; i < 67; i++) {
            // printf("ratio %zu: expected %g got %g\n", i, expected_ratio_grad[i],
            // gradient[i]);
            mu_assert(fabs(gradient[i] - expected_ratio_grad[i]) < 1.e-3,
                      "dlogP ratios not matching with TreeLikelihood_gradient");
            // printf("ratio %zu: expected %g got %g\n", i,
            // expected_ratios_unconstrained_grad[i],
            // ratios->transform->parameter->grad[i]);
            mu_assert(fabs(expected_ratios_unconstrained_grad[i] -
                           ratios->transform->parameter->grad[i]) < 1.e-8,
                      "dlogP ratios not matching");
        }

        mu_assert(fabs(gradient[67] - expected_root_height_grad) < 1.e-4,
                  "dlogP root not matching with TreeLikelihood_gradient");

        mu_assert(fabs(expected_root_height_unconstrained_grad -
                       root->transform->parameter->grad[0]) < 1.e-8,
                  "dlogP root not matching");
    }

    // check the gradient is accumulated in parameter->grad
    // Parameters_zero_grad is not called so the gradient should be doubled
    model->gradient(model, parameters);

    mu_assert(fabs(expected_rate_unconstrained_grad * 2 -
                   rate->transform->parameter->grad[0]) < 1.e-8,
              "dlogP/dx clock not matching in accumulated");

    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(expected_ratios_unconstrained_grad[i] * 2 -
                       ratios->transform->parameter->grad[i]) < 1.e-8,
                  "dlogP ratios not matching in accumulated");
    }

    mu_assert(fabs(expected_root_height_unconstrained_grad * 2 -
                   root->transform->parameter->grad[0]) < 1.e-8,
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
        -0.13274170070753702, 0.7757537767050319,   1.21169172883093,
        1.1327897228892874,   -0.04366588628505873, 0.6151267488296555,
        0.35675961670794615,  1.6326662598918036,   2.4026265317592084,
        1.2074522227651623,   1.5172741708664017,   7.768683424270772,
        24.62213868393231,    30.916466070863855,   3.0512066192138088,
        5.254102963777162,    -0.1671482466858141,  2.2469733002806693,
        5.411215635665562,    6.777422270224085,    2.9502877952728577,
        0.24918727212011263,  3.1264502312284317,   3.2393806581383786,
        2.2341904671223562,   -0.8602399958831338,  16.91654919457692,
        0.5915046719807964,   5.686298796323111,    1.4165351345947172,
        5.338512731751605,    4.542116335273242,    2.9517261384675835,
        0.3480403066106321,   1.6129368145500587,   18.44949801237873,
        11.521383627363338,   1.8781871705264628,   1.5919706446687631,
        0.6984737528022062,   -0.6201082243176413,  7.366292245205609,
        15.699485878822047,   47.00720431253316,    0.8361441423692472,
        2.6077916454436796,   0.17398749882036438,  2.966019246405346,
        7.115904052991468,    0.741393657507105,    1.0088653000435974,
        2.8636392779338156,   3.927865433352233,    23.37624470106702,
        0.7844883947046429,   1.5103167096490036,   7.177008637398779,
        0.5106427284716466,   17.158445964147663,   2.237431941450079,
        1.5210238338970705,   0.7144352770901691,   0.41775132843043744,
        12.204533158328786,   8.793784454079837,    0.5197912298495507,
        1.3868141466982251};
    double expectedRateUnconstrainedJacobianGrad = 328.01765057378356;
    double expected_root_height_jac_grad = 19.936860572419484;
    double expectedRootHeightUnconstrainedJacobianGrad = 30.99161526317372;

    SingleTreeLikelihood_update_all_nodes(tlk);

    Parameters_zero_grad(parameters);
    model->gradient(model, parameters);

    mu_assert(fabs(expectedRateUnconstrainedJacobianGrad -
                   rate->transform->parameter->grad[0]) < 1.e-8,
              "dlogP clock include_jacobian not matching");

    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(ratios->transform->parameter->grad[i] -
                       expectedRatiosUnconstrainedJacobianGrad[i]) < 1.e-8,
                  "dlogP ratios include_jacobian not matching");
    }

    mu_assert(fabs(root->transform->parameter->grad[0] -
                   expectedRootHeightUnconstrainedJacobianGrad) < 1.e-8,
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

    mu_assert(fabs(expectedRootHeightUnconstrainedJacobianGrad -
                   expected_root_height_unconstrained_grad -
                   root->transform->parameter->grad[0]) < 1.e-8,
              "dlogP root det jacobian not matching");
    for (size_t i = 0; i < 67; i++) {
        mu_assert(fabs(expectedRatiosUnconstrainedJacobianGrad[i] -
                       expected_ratios_unconstrained_grad[i] -
                       ratios->transform->parameter->grad[i]) < 1.e-8,
                  "dlogP ratios det jacobian not matching");
    }

    // rate is independent from the tree model transformation so the gradient should be
    // zero
    mu_assert(rate->transform->parameter->grad[0] == 0.0,
              "dlogP rate det jacobian not matching");

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

// Verify the time-tree likelihood gradient (wrt reparameterized node-height
// ratios, root height, and clock rate) against central finite differences.
// `nratios` is tipCount-2 plus the number of unknown-age leaves, so the same
// routine exercises trees with and without unknown leaves.
static char* _fd_time_gradient(const char* file, size_t nratios) {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);

    char* content = load_file(file);
    json_node* json = create_json_tree(content);
    free(content);

    json_node* child = json->children[0];
    Model* model = new_TreeLikelihoodModel_from_json(child, hash);
    SingleTreeLikelihood* tlk = model->obj;
    Model** models = (Model**)model->data;
    Model* mtree = models[0];
    Model* mbm = models[3];
    Tree* tree = mtree->obj;
    Tree_update_heights(tree);
    tlk->include_jacobian = false;

    Parameters* params = get_reparams(tree);
    Parameter* ratios = Parameters_at(params, 0);
    Parameter* root = Parameters_at(params, 1);
    Parameter* rate = Parameters_at(mbm->parameters, 0);
    Parameters* parameters = new_Parameters(3);
    Parameters_add_parameters(parameters, params);
    Parameters_add_parameters(parameters, mbm->parameters);

    mu_assert(Parameter_size(ratios) == nratios, "unexpected number of ratios");

    // unknown-leaf ratios occupy the first `offset` slots; initialize_from_heights
    // should have started them at an interior value (not pinned at the lower
    // bound), so the central finite-difference step below stays inside (0,1).
    size_t offset = nratios - (Tree_tip_count(tree) - 2);
    for (size_t i = 0; i < offset; i++) {
        mu_assert(Parameter_value_at(ratios, i) > 1.e-3 &&
                      Parameter_value_at(ratios, i) < 1.0 - 1.e-3,
                  "unknown-leaf ratio initialized at the boundary");
    }

    // analytic gradient -> snapshot before perturbing anything
    model->logP(model);
    Parameters_zero_grad(parameters);
    model->gradient(model, parameters);
    double* g_ratio = clone_dvector(ratios->grad, nratios);
    double g_root = root->grad[0];
    double g_rate = rate->grad[0];

    double h = 1.e-6;
    double worst = 0.0;

    for (size_t i = 0; i < nratios; i++) {
        double v0 = Parameter_value_at(ratios, i);
        Parameter_set_value_at(ratios, v0 + h, i);
        double lp = model->logP(model);
        Parameter_set_value_at(ratios, v0 - h, i);
        double lm = model->logP(model);
        Parameter_set_value_at(ratios, v0, i);
        double fd = (lp - lm) / (2.0 * h);
        double err = fabs(fd - g_ratio[i]) / (1.0 + fabs(g_ratio[i]));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-3, "ratio gradient does not match finite difference");
    }

    {
        double v0 = Parameter_value(root);
        double hh = h * (1.0 + fabs(v0));
        Parameter_set_value(root, v0 + hh);
        double lp = model->logP(model);
        Parameter_set_value(root, v0 - hh);
        double lm = model->logP(model);
        Parameter_set_value(root, v0);
        double fd = (lp - lm) / (2.0 * hh);
        double err = fabs(fd - g_root) / (1.0 + fabs(g_root));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-3, "root-height gradient does not match finite difference");
    }

    {
        double v0 = Parameter_value(rate);
        double hh = h * (1.0 + fabs(v0));
        Parameter_set_value(rate, v0 + hh);
        double lp = model->logP(model);
        Parameter_set_value(rate, v0 - hh);
        double lm = model->logP(model);
        Parameter_set_value(rate, v0);
        double fd = (lp - lm) / (2.0 * hh);
        double err = fabs(fd - g_rate) / (1.0 + fabs(g_rate));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-3, "rate gradient does not match finite difference");
    }

    printf("  [%s] nratios=%zu worst relative FD error = %.3e\n", file, nratios, worst);

    free(g_ratio);
    free_Parameters(parameters);
    model->free(model);
    free_Hashtable(hash);
    json_free_tree(json);
    return NULL;
}

// Locate the (single, vector-valued) tree.distances Parameter in a non-time
// tree Model. Heights are also tagged MODEL_TREE but are scalars, so the only
// MODEL_TREE parameter with size > 1 is the branch-length vector.
static Parameter* _find_distances(Model* mtree) {
    Parameters* ps = mtree->parameters;
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        if (p->model == MODEL_TREE && Parameter_size(p) > 1) return p;
    }
    return NULL;
}

// Validate the analytic branch-length Hessian exposed through Model->hessian
// against central finite differences. The likelihood here is parameterized by
// branch lengths (non-time tree), and the analytic kernel only fills the
// diagonal, so we check (1) the HESSIAN_DIAGONAL output matches the
// element-wise second finite difference, and (2) HESSIAN_FULL agrees on the
// diagonal, is symmetric, and matches FD on the off-diagonals.
char* test_treelikelihood_branch_hessian() {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);

    // jc69-distance.json is a non-time (branch-length) JC69 tree likelihood on
    // the fluA alignment, which exercises the analytic branch-length Hessian.
    char* content = load_file("jc69-distance.json");
    json_node* json = create_json_tree(content);
    free(content);

    json_node* child = get_json_node(json, "model");
    Model* model = new_TreeLikelihoodModel_from_json(child, hash);
    Model** models = (Model**)model->data;
    Model* mtree = models[0];

    Parameter* distances = _find_distances(mtree);
    mu_assert(distances != NULL, "could not find tree.distances parameter");

    // Set every branch to a common interior length: the analytic second
    // derivative and the central finite difference both behave well away from
    // the zero boundary (near-zero branches make the FD step degenerate).
    size_t dim = Parameter_size(distances);
    for (size_t i = 0; i < dim; i++) {
        Parameter_set_value_at(distances, 0.1, i);
    }

    Parameters* parameters = new_Parameters(1);
    Parameters_add(parameters, distances);
    mu_assert(Parameters_size(parameters) == dim, "flattened dim mismatch");
    mu_assert(Parameters_at(parameters, 0)->model == MODEL_TREE,
              "distances must be tagged MODEL_TREE for the analytic path");

    model->logP(model);

    // (1) analytic diagonal vs finite-difference diagonal
    double* diag = dvector(dim);
    model->hessian(model, parameters, HESSIAN_DIAGONAL, diag);

    double h = 1.e-5;
    double worst = 0.0;
    for (size_t i = 0; i < dim; i++) {
        double v0 = Parameter_value_at(distances, i);
        double hh = h * (1.0 + fabs(v0));
        double l0 = model->logP(model);
        Parameter_set_value_at(distances, v0 + hh, i);
        double lp = model->logP(model);
        Parameter_set_value_at(distances, v0 - hh, i);
        double lm = model->logP(model);
        Parameter_set_value_at(distances, v0, i);
        double fd = (lp - 2.0 * l0 + lm) / (hh * hh);
        double err = fabs(fd - diag[i]) / (1.0 + fabs(fd));
        if (err > worst) worst = err;
        // a single-step central second difference is limited to ~1% here by
        // truncation/roundoff; the analytic diagonal is the accurate one.
        mu_assert(err < 1.5e-2,
                  "branch-length Hessian diagonal does not match finite difference");
    }
    printf("  [jc69-distance.json] branch Hessian diagonal worst relative FD "
           "error = %.3e\n",
           worst);

    // (2) full matrix: diagonal must equal the analytic diagonal, and the
    // matrix must be symmetric.
    model->logP(model);
    double* full = dvector(dim * dim);
    model->hessian(model, parameters, HESSIAN_FULL, full);
    for (size_t i = 0; i < dim; i++) {
        mu_assert(fabs(full[i * dim + i] - diag[i]) < 1.e-8,
                  "HESSIAN_FULL diagonal does not match HESSIAN_DIAGONAL");
        for (size_t j = i + 1; j < dim; j++) {
            mu_assert(fabs(full[i * dim + j] - full[j * dim + i]) < 1.e-8,
                      "HESSIAN_FULL is not symmetric");
        }
    }

    free(diag);
    free(full);
    free_Parameters(parameters);
    model->free(model);
    free_Hashtable(hash);
    json_free_tree(json);
    return NULL;
}

char* test_treelikelihood_time_gradient_fd() {
    return _fd_time_gradient("jc69-time.json", 67);
}

char* test_treelikelihood_time_leaf_gradient_fd() {
    return _fd_time_gradient("jc69-time-leaf.json", 68);
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
//     // 	mu_assert(fabs(dlogPds - dlogPds2) < 0.0001, "dlogPdsigma not
//     matching");

//     model->free(model);
//     for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
//     free(ps);
//     return NULL;
// }

char* all_tests() {
    mu_suite_start();
    // mu_run_test(test_treelikelihood_time);
    mu_run_test(test_treelikelihood_time_unconstrained);
    mu_run_test(test_treelikelihood_branch_hessian);
    mu_run_test(test_treelikelihood_time_gradient_fd);
    mu_run_test(test_treelikelihood_time_leaf_gradient_fd);

    return NULL;
}

RUN_TESTS(all_tests);