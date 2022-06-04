#include "minunit.h"

#include "phyc/filereader.h"
#include "phyc/hashtable.h"
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
    Tree_update_heights(tree);

    tlk->include_jacobian = false;

    double expected_logP = -4777.616349713985;
    double logP = model->logP(model);
    mu_assert(fabs(logP - expected_logP) < 1.e-8, "logP not matching");

    Parameters* ps = new_Parameters(10);
    for (size_t i = 0; i < Parameters_count(bm->rates); i++) {
        Parameters_add(ps, Parameters_at(bm->rates, i));
    }
    model->prepare_gradient(model, ps);
    double expected_rate_grad = 328017.6732813406;
    double dlogP = model->dlogP(model, Parameters_at(ps, 0));
    mu_assert(fabs(dlogP - expected_rate_grad) < 1.e-8, "dlogP clock not matching");

    Parameters* ratios = get_reparams(tree);
    for (size_t i = 0; i < Parameters_count(ratios); i++) {
        Parameters_add(ps, Parameters_at(ratios, i));
    }

    model->prepare_gradient(model, ps);
    SingleTreeLikelihood_update_all_nodes(tlk);
    dlogP = model->dlogP(model, Parameters_at(ps, 0));
    mu_assert(fabs(dlogP - expected_rate_grad) < 1.e-8, "dlogP clock not matching");

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
    for (size_t i = 0; i < Parameters_count(ratios) - 1; i++) {
        dlogP = model->dlogP(model, Parameters_at(ratios, i));
        mu_assert(fabs(dlogP - expected_ratio_grad[i]) < 1.e-8,
                  "dlogP ratios not matching");
    }
    dlogP = model->dlogP(model, Parameters_at(ratios, Parameters_count(ratios) - 1));
    mu_assert(fabs(dlogP - expected_root_height_grad) < 1.e-8,
              "dlogP root not matching");

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

    dlogP = model->dlogP(model, Parameters_at(ps, 0));
    mu_assert(fabs(dlogP - expected_rate_grad) < 1.e-8,
              "dlogP clock include_jacobian not matching");

    for (size_t i = 0; i < Parameters_count(ratios) - 1; i++) {
        dlogP = model->dlogP(model, Parameters_at(ratios, i));
        mu_assert(fabs(dlogP - expected_ratios_jac_grad[i]) < 1.e-8,
                  "dlogP ratios include_jacobian not matching");
    }

    dlogP = model->dlogP(model, Parameters_at(ratios, Parameters_count(ratios) - 1));
    mu_assert(fabs(dlogP - expected_root_height_jac_grad) < 1.e-8,
              "dlogP root include_jacobian not matching");

    free_Parameters(ps);
    model->free(model);
    free_Hashtable(hash);
    json_free_tree(json);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_treelikelihood_time);

    return NULL;
}

RUN_TESTS(all_tests);