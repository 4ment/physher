//
//  vartreelikelihood.c
//  physher
//
//  Created by mathieu on 16/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "vartreelikelihood.h"

#include "distmodel.h"
#include "matrix.h"

VariationalTreeLikelihood * new_VariationalTreeLikelihood(Tree* tree, BranchModel* bm, DistributionModel* distribution){
    VariationalTreeLikelihood* var = malloc(sizeof(VariationalTreeLikelihood));
    var->tree = tree;
    var->bm = bm;
    var->distribution = distribution;
    var->update = true;
    var->gradient = dvector(Tree_tip_count(tree)-1+ Parameters_count(bm->rates));
    var->gradient_distance = dvector(Tree_node_count(tree));
    return var;
}

void _varTreelikelihood_handle_change( Model *self, Model *model, int index ){
    VariationalTreeLikelihood *tlk = self->obj;
    
    if ( model->type == MODEL_TREE ) {
        tlk->update = true;
    }
    else if ( model->type == MODEL_BRANCHMODEL ) {
        tlk->update = true;
    }
    else {
        fprintf(stderr, "%s of type %s\n", model->name, model_type_strings[model->type]);
        error("Unknown change in _varTreelikelihood_handle_change\n");
    }
}

void _varTreeLikelihood_update(VariationalTreeLikelihood* var){
    Tree* tree = var->tree;
    Node** nodes = Tree_nodes(tree);
    Node* root = Tree_root(tree);
    Node* left = root->left;
    Node* right = root->right;
    size_t node_count = Tree_node_count(tree);
    for (size_t i = 0; i < node_count; i++) {
        if(nodes[i] != root && !Node_isroot(Node_parent(nodes[i]))){
//            double blens = Node_time_elapsed(nodes[i]) * var->bm->get(var->bm, nodes[i]);
            double blens = Tree_node_time_elapsed(tree, nodes[i]) * var->bm->get(var->bm, nodes[i]);
            Node_set_distance(nodes[i], blens);
//            printf("blens %f\n", blens);
        }
    }
//    double blens = Node_time_elapsed(left) * var->bm->get(var->bm, left) + Node_time_elapsed(right) * var->bm->get(var->bm, right);
    double blens = Tree_node_time_elapsed(tree, left) * var->bm->get(var->bm, left) + Tree_node_time_elapsed(tree, right) * var->bm->get(var->bm, right);
    Node_set_distance(left, blens);
    Node_set_distance(right, INFINITY);
//    printf("%f %f\n", blens, Node_height(root));
}

double _varTreeLikelihood_logP(Model *self){
    VariationalTreeLikelihood* var = self->obj;
    if (var->update) {
        _varTreeLikelihood_update(var);
        var->update = false;
    }
    return var->distribution->logP(var->distribution);
}

static double _epoch_gradient_addition(const Node* node, const Node* child, const double* lowers, const double* ratios, const double* ratio_gradient) {
    if (Node_isleaf(child)) {
        return 0.0;
    }
    size_t node_id = Node_id(node);
    size_t child_id = Node_id(child);
    size_t node_class_id = Node_class_id(node);
    size_t child_class_id = Node_class_id(child);
    
    // child_id and node_id are in the same epoch
    if (lowers[node_id] == lowers[child_id]) {
        return ratio_gradient[child_class_id] * ratios[child_class_id] / ratios[node_class_id];
    }
    // NOT the same epoch
    else {
        double height = Node_height(node);
        return ratio_gradient[child_class_id] * ratios[child_class_id] / (height - lowers[child_id]) * (height - lowers[node_id])/ratios[node_class_id];
    }
}

static double _height_gradient(Tree* tree, const double* ratios, const double* gradient_height) {
    size_t gradient_length = Tree_tip_count(tree) - 1;
    double* multiplierArray = clone_dvector(ratios, gradient_length);
    multiplierArray[Node_class_id(Tree_root(tree))] = 1.0;
    Node** nodes = Tree_get_nodes(tree, PREORDER);
    size_t node_count = Tree_node_count(tree);
    
    for(size_t i = 1; i < node_count; i++){
        if(!Node_isleaf(nodes[i]))
        multiplierArray[Node_class_id(nodes[i])] *= multiplierArray[Node_class_id(nodes[i]->parent)];
    }
    double sum = 0.0;
    for (int i = 0; i < gradient_length; i++) {
        sum += gradient_height[i] * multiplierArray[i];
    }
    free(multiplierArray);
    return sum;
}

void _update_ratios_gradient(Tree* tree, const double* lowers, const double*ratios, const double* gradient_height, double* gradient){
//double* _update_ratios_gradient(Tree* tree, const double* lowers, const double*ratios, const double* gradient_height){
    Node** nodes = Tree_get_nodes(tree, POSTORDER);
//    double* gradient = dvector(Tree_tip_count(tree)-1);
    for (int i = 0; i < Tree_node_count(tree)-1; i++) {
        if(Node_isleaf(nodes[i])) continue;
        int node_id = Node_id(nodes[i]);
        int node_class_id = Node_class_id(nodes[i]);
        double height = Node_height(nodes[i]);
        gradient[node_class_id] += (height - lowers[node_id])/ratios[node_class_id] * gradient_height[node_class_id];
        gradient[node_class_id] += _epoch_gradient_addition(nodes[i], nodes[i]->left, lowers, ratios, gradient);
        gradient[node_class_id] += _epoch_gradient_addition(nodes[i], nodes[i]->right, lowers, ratios, gradient);
    }
//    return gradient;
}

//static void _varTreeLikelihood_gradient(VariationalTreeLikelihood* var){
static void _varTreeLikelihood_gradient(Model* model){
    Model** subs = model->data;
    VariationalTreeLikelihood* var = model->obj;
    Node** nodes = Tree_get_nodes(var->tree, POSTORDER);
    Node* root = Tree_root(var->tree);
    size_t root_class_id = Node_class_id(root);
    size_t node_count = Tree_node_count(var->tree);
    
    double* gradient_height = dvector(Tree_tip_count(var->tree)-1);
//    double* lowers = dvector(node_count);
    double* lowers = Tree_lowers(var->tree);
    double* ratios = dvector(Tree_tip_count(var->tree)-1);
    double* log_time = dvector(Tree_tip_count(var->tree)-1);
    
    unsigned *map = get_reparam_map(var->tree);
    Parameters* reparams = get_reparams(var->tree);
    
    size_t height_grad_len = Tree_tip_count(var->tree)-1;
    size_t grad_len = height_grad_len;
    if(true){
        grad_len += Parameters_count(var->bm->rates);
    }
    memset(var->gradient, 0, sizeof(double)*grad_len);
    
    for (int i = 0; i < node_count-2; i++) {
        var->gradient_distance[Node_id(nodes[i])] = var->distribution->dlogP(var->distribution, nodes[i]->distance);
    }
    var->gradient_distance[Node_id(root->right)] = var->gradient_distance[Node_id(root->left)];
    
    // strict clock
    if(Parameters_count(var->bm->rates) == 1){
        for (int i = 0; i < node_count-1; i++) {
            var->gradient[height_grad_len] += var->gradient_distance[Node_id(nodes[i])] * Node_time_elapsed(nodes[i]);
        }
    }
    // 1 rate per branch clock
    else{
        for (int i = 0; i < node_count-1; i++) {
            var->gradient[height_grad_len+Node_class_id(nodes[i])] = var->gradient_distance[Node_id(nodes[i])] * Node_time_elapsed(nodes[i]);
        }
    }

    for (int i = 0; i < node_count-1; i++) {
        var->gradient_distance[Node_id(nodes[i])] *= var->bm->get(var->bm, nodes[i]);
    }
    
    for (int i = 0; i < node_count; i++) {
        size_t node_class_id = Node_class_id(nodes[i]);
        size_t node_id = Node_id(nodes[i]);
        
        if(!Node_isleaf(nodes[i])){
            if(nodes[i] != root){
                gradient_height[node_class_id] = -var->gradient_distance[node_id];
            }
            size_t node_left_id = Node_id(nodes[i]->left);
            size_t node_right_id = Node_id(nodes[i]->right);
            gradient_height[node_class_id] += var->gradient_distance[node_left_id];
            gradient_height[node_class_id] += var->gradient_distance[node_right_id];
            
            lowers[node_id] = dmax(lowers[node_left_id], lowers[node_right_id]);
            ratios[node_class_id] = Parameters_value(reparams, map[node_id]);
            log_time[node_class_id] = 1.0 / (Node_height(nodes[i]) - lowers[node_id]);
        }
        else{
            lowers[node_id] = Node_height(nodes[i]);
        }
    }
    //
    log_time[root_class_id] = 0;
    
    _update_ratios_gradient(var->tree, lowers, ratios, gradient_height, var->gradient);
    var->gradient[root_class_id] = _height_gradient(var->tree, ratios, gradient_height);
//    double* temp = _update_ratios_gradient(var->tree, lowers, ratios, gradient_height);
//    temp[root_class_id] = _height_gradient(var->tree, ratios, gradient_height);
    
//    double eps = 0.00001;
//
//    for(size_t i = 0; i < node_count; i++){
//        if(!Node_isleaf(nodes[i])){
////            if(!Node_isroot(nodes[i]->parent) && !Node_isroot(nodes[i])) continue;
//            double r = Parameters_value(reparams, map[Node_id(nodes[i])]);
//
//            Parameters_set_value(reparams, map[Node_id(nodes[i])], r -eps);
//            _varTreeLikelihood_update(var);
//            double mm = var->distribution->logP(var->distribution);
//
//            Parameters_set_value(reparams, map[Node_id(nodes[i])], r +eps);
//            _varTreeLikelihood_update(var);
//            double pp = var->distribution->logP(var->distribution);
//
//            gradient_height[Node_class_id(nodes[i])] = (pp-mm)/(2*eps);
////            var->gradient[Node_class_id(nodes[i])] = (pp-mm)/(2*eps);
//            Parameters_set_value(reparams, map[Node_id(nodes[i])], r);
//        }
//    }
//    _varTreeLikelihood_update(var);
//
//    for(size_t i = 0; i < node_count; i++){
//    if(!Node_isleaf(nodes[i])){
//        printf("%s %e %e %e %f %s\n", nodes[i]->name, gradient_height[Node_class_id(nodes[i])], var->gradient[Node_class_id(nodes[i])], ratios[Node_class_id(nodes[i])], lowers[Node_id(nodes[i])], Node_isroot(nodes[i]) ? "" :Node_isroot(nodes[i]->parent) ? "*": "");
//    }}
//    exit(2);
    // Add gradient of log Jacobian determinant
    _update_ratios_gradient(var->tree, lowers, ratios, log_time, var->gradient);
    var->gradient[root_class_id] += _height_gradient(var->tree, ratios, log_time);
//    double* jac = _update_ratios_gradient(var->tree, lowers, ratios, log_time);
//    jac[root_class_id] += _height_gradient(var->tree, ratios, log_time);

    for (int i = 0; i < node_count-1; i++) {
        if(!Node_isleaf(nodes[i])){
            size_t node_class_id = Node_class_id(nodes[i]);
            var->gradient[node_class_id] -= 1.0 / ratios[node_class_id];
        }
    }
//    for (int i = 0; i < node_count-1; i++) {
//        if(!Node_isleaf(nodes[i])){
//            size_t node_class_id = Node_class_id(nodes[i]);
//            var->gradient[node_class_id] = temp[node_class_id] + jac[node_class_id] - 1.0 / ratios[node_class_id];
//        }
//    }
//    var->gradient[root_class_id] = temp[root_class_id] + jac[root_class_id];
//    free(temp);
//    free(jac);
    
//    free(lowers);
    free(ratios);
    free(log_time);
    free(gradient_height);
}

double _varTreeLikelihood_dlogP(Model *self, const Parameter* p){
    VariationalTreeLikelihood* var = self->obj;
    if (var->update) {
        _varTreeLikelihood_update(var);
        _varTreeLikelihood_gradient(self);
        var->update = false;
    }
    if (p->id < 0) {
        int realIndex = -p->id-1;
        Node* node = Tree_node(var->tree, realIndex);
        return var->gradient[Node_class_id(node)];
    }
    else if(p->model == MODEL_BRANCHMODEL){
        if(Parameters_count(var->bm->rates) == 1){
            return var->gradient[Tree_tip_count(var->tree)-1];
        }
        else{
            Node* node = Tree_node(var->tree, p->id);
            return var->gradient[Tree_tip_count(var->tree)-1+Node_class_id(node)];
        }
    }
    return 0;
}

static void _varTreeLikelihood_model_free( Model *self ){
    if(self->ref_count == 1){
        //printf("Free variationalTreelikelihood model %s\n", self->name);
        VariationalTreeLikelihood* var = (VariationalTreeLikelihood*)self->obj;
        Model** list = (Model**)self->data;
        for(int i = 0; i < 3; i++){
            list[i]->free(list[i]);
        }
        free(var->gradient);
        free(var->gradient_distance);
        free(list);
        free(var);
        free_Model(self);
    }
    else{
        self->ref_count--;
    }
}

static Model* _varTreeLikelihood_model_clone(Model* self, Hashtable* hash){
    error("_varTreeLikelihood_model_clone not implemented yet");
    Model** list = (Model**)self->data;
    Model* mtree = list[0];
    Model* mbm = list[1];
    Model* mdist = list[2];
    Model *mtreeclone = NULL;
    Model* mbmclone = NULL;
    Model* mdistclone = NULL;
    
    if(Hashtable_exists(hash, mtree->name)){
        mtreeclone = Hashtable_get(hash, mtree->name);
        mtreeclone->ref_count++;
    }
    else{
        mtreeclone = mtree->clone(mtree, hash);
        Hashtable_add(hash, mtreeclone->name, mtreeclone);
    }
    
    if(Hashtable_exists(hash, mbm->name)){
        mbmclone = Hashtable_get(hash, mbm->name);
        mbmclone->ref_count++;
    }
    else{
        mbmclone = mbm->clone(mbm, hash);
        Hashtable_add(hash, mbmclone->name, mbmclone);
    }
        
    if(Hashtable_exists(hash, mdist->name)){
        mdistclone = Hashtable_get(hash, mdist->name);
        mdistclone->ref_count++;
    }
    else{
        mdistclone = mtree->clone(mdist, hash);
        Hashtable_add(hash, mdistclone->name, mdistclone);
    }
    
    VariationalTreeLikelihood* var = (VariationalTreeLikelihood*)self->obj;
    VariationalTreeLikelihood* clonevar = NULL;
    Model* clone = NULL;
    

//    clonevar = clone_SingleTreeLikelihood_with(tlk, (Tree*)mtreeclone->obj, (SiteModel*)msmclone->obj, tlk->sp, NULL);
//    clone = new_TreeLikelihoodModel(self->name, clonevar, mtreeclone, mbmclone, mdistclone);
    
    Hashtable_add(hash, clone->name, clone);
    mtreeclone->free(mtreeclone);
    mbmclone->free(mbmclone);
    mdistclone->free(mdistclone);
    clone->store = self->store;
    clone->restore = self->restore;
    clone->storedLogP = self->storedLogP;
    clone->lp = self->lp;
    clone->full_logP = self->full_logP;
    return clone;
}

Model * new_VariationalTreeLikelihoodModel( const char* name, VariationalTreeLikelihood *tlk,  Model *tree, Model *bm, Model* distribution ){
    Model *model = new_Model(MODEL_VARIATIONAL_TREELIKELIHOOD, name, tlk);

    tree->listeners->add( tree->listeners, model );
    bm->listeners->add( bm->listeners, model );
    
    model->logP = _varTreeLikelihood_logP;
    model->full_logP = _varTreeLikelihood_logP;
    model->dlogP = _varTreeLikelihood_dlogP;
    model->d2logP = NULL;
    model->ddlogP = NULL;
    model->handle_restore = NULL;
    model->store = NULL;
    model->restore = NULL;
    
    model->update = _varTreelikelihood_handle_change;
    model->free = _varTreeLikelihood_model_free;
    model->clone = _varTreeLikelihood_model_clone;
    Model** list = (Model**)malloc(sizeof(Model*)*3);
    model->data = list;
    list[0] = tree;
    list[1] = bm;
    list[2] = distribution;
    tree->ref_count++;
    bm->ref_count++;
    distribution->ref_count++;
    return model;
}

Model* safe_get_model_reference(Hashtable* hash, char* ref){
    if(strlen(ref) <= 1){
        fprintf(stderr, "invlid ref: %s\n", ref);
        exit(1);
    }
    if(ref[0] != '&'){
        fprintf(stderr, "ref does not start with &: %s\n", ref);
        exit(1);
    }
    Model* m = Hashtable_get(hash, ref+1);
    if(m == NULL){
        fprintf(stderr, "ref %s does not exist\n", ref+1);
        exit(1);
    }
    return m;
}

Model * new_VariationalTreeLikelihoodModel_from_json(json_node*node, Hashtable*hash){
    char* allowed[] = {
        "branchmodel",
        "distribution",
        "tree"
    };
    json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
    
    json_node* bm_node = get_json_node(node, "branchmodel");
    json_node* dist_node = get_json_node(node, "distribution");
    json_node* tree_node = get_json_node(node, "tree");
    
    Model* mtree = NULL;
    Model* mbm = NULL;
    Model* mdist = NULL;
    
    if (tree_node->node_type == MJSON_STRING) {
        mtree = safe_get_model_reference(hash, tree_node->value);
        mtree->ref_count++;
    }
    else{
        char* id = get_json_node_value_string(tree_node, "id");
        mtree = new_TreeModel_from_json(tree_node, hash);
        Hashtable_add(hash, id, mtree);
    }
    
    if (bm_node->node_type == MJSON_STRING) {
        mbm = safe_get_model_reference(hash, bm_node->value);
        mbm->ref_count++;
    }
    else{
        char* id = get_json_node_value_string(bm_node, "id");
        mbm = new_BranchModel_from_json(bm_node, hash);
        Hashtable_add(hash, id, mbm);
    }
    
    
    if (dist_node->node_type == MJSON_STRING) {
        mdist = safe_get_model_reference(hash, dist_node->value);
        mdist->ref_count++;
    }
    else{
        char* id = get_json_node_value_string(dist_node, "id");
        mdist = new_DistributionModel_from_json(dist_node, hash);
        Hashtable_add(hash, id, mdist);
    }
    
    VariationalTreeLikelihood* tlk = new_VariationalTreeLikelihood(mtree->obj, mbm->obj, mdist->obj);
    char* id = get_json_node_value_string(node, "id");
    Model* model = new_VariationalTreeLikelihoodModel(id, tlk, mtree, mbm, mdist);
    mtree->free(mtree);
    mbm->free(mbm);
    mdist->free(mdist);

    return model;
}
