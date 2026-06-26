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

//TODO: check if it is used and how
static void _jacobian_model_gradient(Model* self, Parameters* parameters) {
    Parameters* ps = self->obj;
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        Transform* t = p->transform;
        t->gradient_log_det_jacobian(t);
    }
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

static void _jacobian_model_store(Model* self){
    if(!self->stored){
        self->storedLogP = self->lp;
        Parameters* parameters = self->obj;
        Parameters_store(parameters);
        self->stored = true;
    }
}

static void _jacobian_model_restore(Model* self){
    if(self->stored){
        self->lp = self->storedLogP;
        Parameters* parameters = self->obj;
        Parameters_restore(parameters);
        self->stored = false;
    }
}

static void _jacobian_model_accept(Model* self){
    if(self->stored){
        Parameters* parameters = self->obj;
        Parameters_accept(parameters);
        self->stored = false;
    }
}

Model* new_JacobianTransformModel(const char* id, Parameters* parameters) {
    Model* model = new_Model(MODEL_JACOBIAN_TRANSFORM, id, parameters);
    model->logP = _jacobian_model_logP;
    model->full_logP = _jacobian_model_logP;
    model->gradient = _jacobian_model_gradient;
    model->store = _jacobian_model_store;
    model->restore = _jacobian_model_restore;
    model->accept = _jacobian_model_accept;
    model->free = _jacobian_model_free;
    return model;
}

Model* new_JacobianTransformModel_from_json(json_node* node, Hashtable* hash) {
    static const json_field schema[] = {
        {"parameters", JSON_REQUIRED, JSON_OBJECT_OR_STRING},
    };
    json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
    json_node* parametersNode = get_json_node(node, "parameters");
    Parameters* parameters = new_Parameters(1);
    grab_parameters(parametersNode, hash, parameters);
    const char* id = get_json_node_value_string(node, "id");
    Model* model = new_JacobianTransformModel(id, parameters);
    Hashtable_add(hash, id, model);
    return model;
}