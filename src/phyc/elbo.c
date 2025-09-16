
#include "elbo.h"

#include <stdlib.h>
#include <strings.h>

#include "compoundmodel.h"
#include "distmodelfactory.h"
#include "matrix.h"
#include "treelikelihood.h"
#include "bound.h"

double _elbo_logP(Bound* elbo) {
    size_t infCount = 0;
    double sumLogP = 0;
    double sumLogQ = 0;
    size_t paramCount = Parameters_count(elbo->parameters);
    Model** distributions = &elbo->variational;
    size_t distCount = 1;
    if (elbo->variational->type == MODEL_COMPOUND) {
        CompoundModel* cm = elbo->variational->obj;
        distributions = cm->models;
        distCount = cm->count;
    }

    double entropy = 0;
    // analytical entropy
    if (elbo->entropy) {
        for (size_t i = 0; i < distCount; i++) {
            DistributionModel* dist = distributions[i]->obj;
            if(dist->entropy != NULL){
                entropy += dist->entropy(dist);
            }
        }
    }

    for (size_t i = 0; i < elbo->samples; i++) {
        elbo->variational->rsample(elbo->variational);
        double logP = elbo->joint->logP(elbo->joint);

        if (!isinf(logP) && !isnan(logP)) {
            sumLogP += logP;
            if (!elbo->entropy) {
                sumLogQ += elbo->variational->logP(elbo->variational);
            }

        } else {
            infCount++;
            i--;
        }
        if (infCount == 10) {
            return logP;
        }
    }
    if (!elbo->entropy) {
        entropy = -sumLogQ / (elbo->samples - infCount);
    }
    return sumLogP / (elbo->samples - infCount) + entropy;
}

double _elbo_log_multi(Bound* elbo) {
    size_t infCount = 0;
    size_t paramCount = Parameters_count(elbo->parameters);
    Model** distributions = &elbo->variational;
    size_t distCount = 1;
    if (elbo->variational->type == MODEL_COMPOUND) {
        CompoundModel* cm = elbo->variational->obj;
        distributions = cm->models;
        distCount = cm->count;
    }

    double* elbos = dvector(elbo->kSamples);
    double elboValue = 0;
    for (size_t i = 0; i < elbo->samples; i++) {
        for (size_t k = 0; k < elbo->kSamples; k++) {
            elbo->variational->rsample(elbo->variational);
            double logQ = elbo->variational->logP(elbo->variational);
            double logP = elbo->joint->logP(elbo->joint);

            if (!isinf(logP) && !isnan(logP)) {
                elbos[k] = logP - logQ;
            } else {
                infCount++;
                k--;
            }
            if (infCount == 10) {
                free(elbos);
                return NAN;
            }
        }
        double maxLogELBO = dmax_vector(elbos, elbo->kSamples);
        double mean = 0;
        for (size_t k = 0; k < elbo->kSamples - infCount; k++) {
            mean += exp(elbos[k] - maxLogELBO);
        }
        mean /= elbo->kSamples - infCount;
        elboValue += log(mean) + maxLogELBO;
    }
    free(elbos);
    return elboValue / elbo->samples;
}

double _elbo_gradient(Bound* self, Parameters* parameters) {
    size_t varTotalSize = Parameters_size(parameters);
    size_t varParameterCount = Parameters_count(parameters);
    Model** distributions = &self->variational;
    size_t distCount = 1;
    if (self->variational->type == MODEL_COMPOUND) {
        CompoundModel* cm = self->variational->obj;
        distributions = cm->models;
        distCount = cm->count;
    }

    for (size_t i = 0; i < self->samples; i++) {
        Parameters_zero_grad(self->parameters);

        // sample from variational distribution using reparameterization trick
        self->variational->rsample(self->variational);
        double logQ = self->variational->logP(self->variational);
        double logP = self->joint->logP(self->joint);

        // calculate gradient of joint wrt phylogenetic parameters z
        // z are unconstrained or match the support of the variational distributions
        self->joint->gradient(self->joint, self->parameters);
        // dm->x->grad contains dP/dz

        for (size_t j = 0; j < distCount; j++) {
            DistributionModel* dm = distributions[j]->obj;

            // calculate gradient of joint wrt variational parameters φ using the chain rule arising from the reparameterization trick
            // use dm->x->grad to calculate dP/dφ = dP/dz dz/dφ
            dm->rgradient(dm);
            // dm->parameters->grad contains dP/dφ

            if (!self->entropy || dm->gradient_entropy == NULL) {
                size_t varParameterCount = Parameters_count(dm->parameters);
                Parameters* varParameters = new_Parameters(varParameterCount);
                double* gradP = dvector(Parameters_size(dm->parameters));

                Parameters_zero_grad(dm->x);

                size_t index = 0;
                for(size_t k = 0; k < varParameterCount; k++){
                    Parameter* p = Parameters_at(dm->parameters, k);
                    Parameter *px = Parameters_depends(parameters, p);
                    // check the parameter is not fixed
                    if(px != NULL){
                        memcpy(gradP+index, px->grad, Parameter_size(px)*sizeof(double));
                        index += Parameter_size(px);
                        Parameter_zero_grad(px);
                        Parameters_add(varParameters, px);
                    }
                }
                // Parameters_add_parameters(varParameters, dm->x);
                
                // calculate gradient of variational distribution wrt phylogenetic parameters z
                // z are unconstrained or match the support of the variational distributions
                dm->gradient2(dm, dm->x);
                // dm->x->grad contains dQ/dφ

                // calculate gradient of variational distribution wrt variational parameters φ using the chain rule arising from the reparameterization trick
                // use dm->x->grad to calculate dQ/dφ = dQ/dz dz/dφ
                dm->rgradient(dm);
                // dm->parameters->grad contains dP/dφ

                index = 0;
                for(size_t k = 0; k < Parameters_count(varParameters); k++){
                    Parameter* px = Parameters_at(varParameters, k);
                    for(size_t m = 0; m < Parameter_size(px); m++){
                        px->grad[m] = gradP[index] - px->grad[m];
                        index++;
                    }
                }
                free(gradP);
                free_Parameters(varParameters);
            }
        }
    }

    for (size_t i = 0; i < Parameters_count(parameters); i++) {
        Parameter* p = Parameters_at(parameters, i);
        for (size_t j = 0; j < Parameter_size(p); j++) {
            p->grad[j] /= self->samples;
        }
    }
    
    // add gradient of entropy
    for (size_t j = 0; j < distCount; j++) {
        DistributionModel* dm = distributions[j]->obj;
        if (self->entropy && dm->gradient_entropy != NULL) {
            dm->gradient_entropy(dm, parameters);
        }
    }
    // TODO
    return 0;
}

Model* new_ELBO_from_json(json_node* node, Hashtable* hash) {
    Model* model = new_AbstractBoundModel_from_json(node, hash);
    Bound* bound = model->obj;
    bound->logP = _elbo_logP;
    bound->gradient = _elbo_gradient;
    return model;
}