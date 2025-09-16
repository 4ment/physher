#include "jacobiantransform.h"

#include <stdlib.h>

static double _jacobian_model_logP(Model* self) {
    Parameters* parameters = self->obj;
    double logP = 0;
    for (size_t i = 0; i < Parameters_count(parameters); i++) {
        Parameter* p = Parameters_at(parameters, i);
        Transform* t = p->transform;
        logP += t->log_det_jacobian(t);
    }
    self->lp = logP;
    return logP;
}

static double _jacobian_model_gradient(Model* self, const Parameters* parameters) {
    Parameters* ps = self->obj;
    double sumGrad = 0;
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        Transform* t = p->transform;
        sumGrad += t->gradient_log_det_jacobian(t);
    }
    return sumGrad;
}

static void _jacobian_model_free(Model* self) {
    if(self->ref_count == 1){
        Parameters* parameters = self->obj;
        free_Parameters(parameters);
        free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

Model* new_JacobianTransformModel(const char* id, Parameters* parameters) {
    Model* model = new_Model(MODEL_JACOBIAN_TRANSFORM, id, parameters);
    model->logP = _jacobian_model_logP;
    model->full_logP = _jacobian_model_logP;
    model->gradient = _jacobian_model_gradient;
    model->free = _jacobian_model_free;
    return model;
}

Model* new_JacobianTransformModel_from_json(json_node* node, Hashtable* hash) {
    char* allowed[] = {"parameters"};
    json_check_allowed(node, allowed, sizeof(allowed) / sizeof(allowed[0]));
    json_node* parametersNode = get_json_node(node, "parameters");
    Parameters* parameters = new_Parameters(1);
    grab_parameters(parametersNode, hash, parameters);
    const char* id = get_json_node_value_string(node, "id");
    Model* model = new_JacobianTransformModel(id, parameters);
    Hashtable_add(hash, id, model);
    return model;
}