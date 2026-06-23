
#include "bound.h"

#include <stdlib.h>
#include <strings.h>

#include "compoundmodel.h"
#include "distmodelfactory.h"
#include "matrix.h"
#include "treelikelihood.h"


void free_bound_t(Bound* obj) {
    obj->joint->free(obj->joint);
    obj->variational->free(obj->variational);
    free_Parameters(obj->parameters);
    free(obj);
}

Bound* new_Bound(Model* joint, Model* variational, size_t samples, bool entropy) {
    Bound* obj = malloc(sizeof(Bound));
    obj->logP = NULL;
    obj->gradient = NULL;
    obj->joint = joint;
    obj->variational = variational;
    // obj->joint->ref_count++;
    // obj->variational->ref_count++;
    obj->samples = samples;
    obj->entropy = entropy;
    obj->parameters = new_Parameters(1);

    Model** distributions = &obj->variational;
    size_t distCount = 1;
    if (obj->variational->type == MODEL_COMPOUND) {
        CompoundModel* cm = obj->variational->obj;
        distributions = cm->models;
        distCount = cm->count;
    }

    for (size_t j = 0; j < distCount; j++) {
        DistributionModel* dm = distributions[j]->obj;
        Parameters* ps = dm->x;
        // there is should be only one parameter
        Parameters_add_parameters(obj->parameters, ps);
        // TODO: call prepare gradient of joint here
    }

    return obj;
}

static double _BoundModel_logP(Model* self) {
    Bound* bound = self->obj;
    self->lp = bound->logP(bound);
    return self->lp;
}
static double _BoundModel_full_logP(Model* self) {
    Bound* bound = self->obj;
    self->lp = bound->logP(bound);
    return self->lp;
}

static void _BoundModel_gradient(Model* self, Parameters* parameters) {
    Bound* bound = self->obj;
    bound->gradient(bound, parameters);
}

static void _BoundModel_handle_change(Model* self, Model* model, Parameter* parameter,
                                      int index) {
    // The bound recomputes its objective eagerly (no dirty flag), so there is
    // nothing to invalidate here; re-fire so any model listening to the bound is
    // notified that one of its variational/joint dependencies changed.
    self->listeners->fire(self->listeners, self, parameter, index);
}

static void _BoundModel_handle_restore(Model* self, Model* model, int index) {
    self->listeners->fire_restore(self->listeners, self, index);
}

static void _BoundModel_free(Model* self) {
    if(self->ref_count == 1){
        Bound* bound = self->obj;
        free_bound_t(bound);
        free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

Model* new_BoundModel(const char* name, Bound* bound) {
    Model* model = new_Model(MODEL_BOUNDMODEL, name, bound);
    model->logP = _BoundModel_logP;
    model->full_logP = _BoundModel_full_logP;
    model->gradient = _BoundModel_gradient;
    model->free = _BoundModel_free;
    model->clone = NULL;
    model->store = NULL;
    model->restore = NULL;
    model->sample = NULL;
    model->samplable = false;
    model->update = _BoundModel_handle_change;
    model->handle_restore = _BoundModel_handle_restore;

    // Listen to the submodels so changes to the variational parameters (and the
    // joint's parameters) propagate up to anything listening to the bound. The
    // variational distribution already listens to its own x parameters, so we
    // register on the submodels rather than the parameters to avoid double firing.
    bound->joint->listeners->add(bound->joint->listeners, model);
    bound->variational->listeners->add(bound->variational->listeners, model);

    // self->parameters must hold every parameter the bound depends on; bound->parameters
    // are the variational x parameters collected in new_Bound.
    Parameters_add_parameters_recursively(model->parameters, bound->parameters);
    return model;
}

Model* new_AbstractBoundModel_from_json(json_node* node, Hashtable* hash) {
    char* allowed[] = {"bound", "entropy",  // use entropy or not
                       "joint", "samples", "variational"};
    json_check_allowed(node, allowed, sizeof(allowed) / sizeof(allowed[0]));

    const char* id = get_json_node_value_string(node, "id");

    bool entropy = get_json_node_value_bool(node, "entropy", true);
    size_t samples = 100;
    size_t kSamples = 1;
    json_node* samplesNode = get_json_node(node, "samples");
    if (samplesNode != NULL && samplesNode->node_type == MJSON_ARRAY) {
        int temp[2];
        get_json_node_value_array_int(node, "samples", temp);
        samples = temp[0];
        kSamples = temp[1];
    } else {
        samples = get_json_node_value_size_t(node, "samples", 100);
    }

    json_node* jointNode = get_json_node(node, "joint");
    Model* joint = NULL;
    char* ref = (char*)jointNode->value;
    if (safe_is_reference(ref, id)) {
        joint = safe_get_reference_parameter(ref, hash, id);
        joint->ref_count++;
    } else {
        char* type = get_json_node_value_string(jointNode, "type");
        if (strcasecmp(type, "treelikelihood")) {
            joint = new_TreeLikelihoodModel_from_json(jointNode, hash);
        } else if (strcasecmp(type, "compound")) {
            joint = new_CompoundModel_from_json(jointNode, hash);
        } else if (strcasecmp(type, "distribution")) {
            joint = new_DistributionModel_from_json(jointNode, hash);
        }
    }

    json_node* variationalNode = get_json_node(node, "variational");
    Model* variational = NULL;
    ref = (char*)variationalNode->value;
    if (safe_is_reference(ref, id)) {
        variational = safe_get_reference_parameter(ref, hash, id);
        variational->ref_count++;
    } else {
        char* type = get_json_node_value_string(variationalNode, "type");
        if (strcasecmp(type, "compound")) {
            variational = new_CompoundModel_from_json(variationalNode, hash);
        } else if (strcasecmp(type, "distribution")) {
            variational = new_DistributionModel_from_json(variationalNode, hash);
        }
    }

    Bound* bound = new_Bound(joint, variational, samples, entropy);
    bound->kSamples = kSamples;
    Model* model = new_BoundModel(id, bound);
    return model;
}