//
//  distweibull.c
//  physher
//
//  Created by mathieu on 20/6/26.
//  Copyright © 2026 Mathieu Fourment. All rights reserved.
//

#include "distweibull.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <strings.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"
#include "transforms.h"

// Weibull distribution parameterized with scale (a) and shape (b).
// GSL uses the same scale/shape parameterization:
//   f(x; a, b) = (b/a) (x/a)^(b-1) exp(-(x/a)^b)
//   gsl_ran_weibull_pdf(x, a, b) with a = scale, b = shape

double DistributionModel_weibull_logP(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    const double* scale = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* shape = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior)
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->lp += log(gsl_ran_weibull_pdf(values[i] - dm->shift, scale[0], shape[0]));
            }
        }
    }
    // multiple parameter and len(x) >= 1 and \sum_i len(x_i) == len(parameter)
    else{
        size_t index = 0;
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->lp += log(gsl_ran_weibull_pdf(values[i] - dm->shift, scale[index], shape[index]));
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

// d/dscale log f = b((x/a)^b - 1)/a
// d/dshape log f = 1/b + log(x/a)(1 - (x/a)^b)
// d/dx     log f = (b-1)/x - b(x/a)^b/x
void DistributionModel_weibull_gradient(DistributionModel* dm, Parameters* parameters){
    Parameter* scale = Parameters_at(dm->parameters, 0);
    Parameter* shape = Parameters_at(dm->parameters, 1);

    const double* scaleValues = Parameter_values(scale);
    const double* shapeValues = Parameter_values(shape);
    size_t mask = -(Parameter_size(scale) != 1);

    Parameter* scalex = Parameters_depends(parameters, scale);
    Parameter* shapex = Parameters_depends(parameters, shape);

    if(scalex != NULL){
        // Accumulate the gradient wrt the constrained scale in dm->tempp across
        // all x. A shared scalar scale (mask == 0) folds every element onto
        // tempp[0], so backward sees the full sum, not just the last element.
        memset(dm->tempp, 0, Parameter_size(scale) * sizeof(double));
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                size_t idx = index & mask;
                double scaleValue = scaleValues[idx];
                double shapeValue = shapeValues[idx];
                double z = (xValues[j] - dm->shift) / scaleValue;
                double zb = pow(z, shapeValue);
                double g = shapeValue * (zb - 1.0) / scaleValue;
                dm->tempp[idx] += g;
                scale->grad[idx] += g;
                index++;
            }
        }
        if(scale != scalex){
            scale->transform->backward(scale->transform, dm->tempp);
        }
    }

    if(shapex != NULL){
        memset(dm->tempp, 0, Parameter_size(shape) * sizeof(double));
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                size_t idx = index & mask;
                double scaleValue = scaleValues[idx];
                double shapeValue = shapeValues[idx];
                double z = (xValues[j] - dm->shift) / scaleValue;
                double zb = pow(z, shapeValue);
                double g = 1.0 / shapeValue + log(z) * (1.0 - zb);
                dm->tempp[idx] += g;
                shape->grad[idx] += g;
                index++;
            }
        }
        if(shape != shapex){
            shape->transform->backward(shape->transform, dm->tempp);
        }
    }

    size_t index = 0;
    for(size_t k = 0; k < Parameters_count(dm->x); k++){
        Parameter* x = Parameters_at(dm->x, k);
        Parameter* xx = Parameters_depends(parameters, x);
        size_t sizeX = Parameter_size(x);
        if (xx != NULL) {
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < sizeX; j++){
                size_t idx = index & mask;
                double scaleValue = scaleValues[idx];
                double shapeValue = shapeValues[idx];
                double xs = xValues[j] - dm->shift;
                double z = xs / scaleValue;
                double zb = pow(z, shapeValue);
                // x_k owns its own grad/transform: index into them locally (j)
                // while index & mask tracks the global position for the params.
                dm->tempp[j] = (shapeValue - 1.0) / xs - shapeValue * zb / xs;
                x->grad[j] += dm->tempp[j];
                index++;
            }
            if(x != xx){
                x->transform->backward(x->transform, dm->tempp);
            }
        }
        else{
            index += sizeX;
        }
    }
}

static void DistributionModel_weibull_sample(DistributionModel* dm){
    const double* scale = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* shape = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior)
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_weibull(dm->rng, scale[0], shape[0]) + dm->shift;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
    // multiple parameter and len(x) >= 1 and \sum_i len(x_i) == len(parameter)
    else{
        size_t index = 0;
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_weibull(dm->rng, scale[index], shape[index]) + dm->shift;
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

// Reparameterization trick: eta ~ Uniform(0, 1), x = scale * (-log(1-eta))^(1/shape)
static void DistributionModel_weibull_rsample(DistributionModel* dm){
    const double* scale = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* shape = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);
    size_t index = 0;
    for(size_t j = 0; j < dimX; j++){
        Parameter* x = Parameters_at(dm->x, j);
        size_t sizeX = Parameter_size(x);
        double* temp = dvector(sizeX);
        for (size_t i = 0; i < sizeX; i++) {
            size_t idx = index & -(Parameter_size(Parameters_at(dm->parameters, 0)) != 1);
            dm->tempx[index] = gsl_ran_flat(dm->rng, 0.0, 1.0);  // eta
            temp[i] = scale[idx] * pow(-log(1.0 - dm->tempx[index]), 1.0 / shape[idx]) + dm->shift;
            index++;
        }
        Parameter_set_values(x, temp);
        free(temp);
    }
}

// Gradient of the log PDF wrt the parameters where x was sampled using the
// reparameterization trick (see DistributionModel_weibull_rsample).
//   x = scale * w^(1/shape),  w = -log(1 - eta)
//   dx/dscale = x/scale
//   dx/dshape = -x*log(x/scale)/shape
static void DistributionModel_weibull_rgradient(DistributionModel* dm){
    Parameter* scale = Parameters_at(dm->parameters, 0);
    Parameter* shape = Parameters_at(dm->parameters, 1);
    const double* scaleValues = Parameter_values(scale);
    const double* shapeValues = Parameter_values(shape);
    size_t dimX = Parameters_count(dm->x);
    size_t index = 0;
    for (size_t j = 0; j < dimX; j++) {
        Parameter* x = Parameters_at(dm->x, j);
        size_t dim = Parameter_size(x);
        const double* xValues = Parameter_values(x);
        size_t base = index;
        // dL/dscale = dL/dx * dx/dscale,  dx/dscale = (x - shift)/scale = z
        for (size_t i = 0; i < dim; i++) {
            double z = (xValues[i] - dm->shift) / scaleValues[index];
            dm->tempp[i] = x->grad[i] * z;
            scale->grad[index] += dm->tempp[i];
            index++;
        }
        // Both scale and shape are positive (constrained): push the reparameterized
        // gradient to the unconstrained leaf through the transform when present.
        if (scale->transform != NULL) {
            scale->transform->backward(scale->transform, dm->tempp);
        }
        // dL/dshape = dL/dx * dx/dshape,  dx/dshape = -(x - shift) log(z)/shape
        index = base;
        for (size_t i = 0; i < dim; i++) {
            double xs = xValues[i] - dm->shift;
            double z = xs / scaleValues[index];
            dm->tempp[i] = x->grad[i] * (-xs * log(z) / shapeValues[index]);
            shape->grad[index] += dm->tempp[i];
            index++;
        }
        if (shape->transform != NULL) {
            shape->transform->backward(shape->transform, dm->tempp);
        }
    }
}

// Entropy: H(X) = euler*(1 - 1/shape) + log(scale/shape) + 1
static double DistributionModel_weibull_entropy(DistributionModel* dm){
    const double* scaleValues = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* shapeValues = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimP = Parameter_size(Parameters_at(dm->parameters, 0));
    double entropy = 0;
    for (size_t i = 0; i < dimP; i++) {
        entropy += M_EULER * (1.0 - 1.0 / shapeValues[i]) + log(scaleValues[i] / shapeValues[i]) + 1.0;
    }
    return entropy;
}

static void DistributionModel_weibull_entropy_gradient(DistributionModel* dm, const Parameters* parameters){
    Parameter* scale = Parameters_at(dm->parameters, 0);
    Parameter* shape = Parameters_at(dm->parameters, 1);
    const double* scaleValues = Parameter_values(scale);
    const double* shapeValues = Parameter_values(shape);

    Parameter* scalex = Parameters_depends(parameters, scale);
    if (scalex != NULL) {
        size_t dimP = Parameter_size(scale);
        for (size_t i = 0; i < dimP; i++) {
            dm->tempp[i] = 1.0 / scaleValues[i];
            scale->grad[i] += dm->tempp[i];
        }
        if (scale != scalex) {
            scale->transform->backward(scale->transform, dm->tempp);
        }
    }

    Parameter* shapex = Parameters_depends(parameters, shape);
    if (shapex != NULL) {
        size_t dimP = Parameter_size(shape);
        for (size_t i = 0; i < dimP; i++) {
            dm->tempp[i] = M_EULER / (shapeValues[i] * shapeValues[i]) - 1.0 / shapeValues[i];
            shape->grad[i] += dm->tempp[i];
        }
        if (shape != shapex) {
            shape->transform->backward(shape->transform, dm->tempp);
        }
    }
}

DistributionModel* new_WeibullDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_WEIBULL;
    dm->parameterization = DISTRIBUTION_WEIBULL_SCALE_SHAPE;
    dm->logP = DistributionModel_weibull_logP;
    dm->gradient = DistributionModel_weibull_gradient;
    dm->rgradient = DistributionModel_weibull_rgradient;
    dm->sample = DistributionModel_weibull_sample;
    dm->rsample = DistributionModel_weibull_rsample;
    dm->entropy = DistributionModel_weibull_entropy;
    dm->gradient_entropy = DistributionModel_weibull_entropy_gradient;
    dm->shift = 0;
    dm->support[0] = 0;
    dm->support[1] = INFINITY;
    return dm;
}

Model* new_WeibullDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");

    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));

    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    Parameter* scale = NULL;
    Parameter* shape = NULL;

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* scaleValues = malloc(sizeof(double)*paramCount);
        double* shapeValues = malloc(sizeof(double)*paramCount);

        for (size_t i = 0; i < paramCount; i++) {
            const double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            // Method-of-moments-ish initialization: a moderate shape and a
            // scale matching the sample mean (mean = scale*Gamma(1+1/shape)).
            shapeValues[i] = (v > 0.0 ? m / sqrt(v) : 1.0);
            scaleValues[i] = (m > 0.0 ? m : 1.0);
            free_Vector(samples[i]);
        }

        scale = new_Parameter2("scale", scaleValues, paramCount, new_Constraint(0, INFINITY));
        shape = new_Parameter2("shape", shapeValues, paramCount, new_Constraint(0, INFINITY));

        free(scaleValues);
        free(shapeValues);
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        double* scaleValues = malloc(sizeof(double)*paramCount);
        double* shapeValues = malloc(sizeof(double)*paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            scaleValues[i] = 1.0;
            shapeValues[i] = 1.0;
        }
        scale = new_Parameter2("scale", scaleValues, paramCount, new_Constraint(0, INFINITY));
        shape = new_Parameter2("shape", shapeValues, paramCount, new_Constraint(0, INFINITY));

        free(scaleValues);
        free(shapeValues);
    }
    else{
        json_node* parametersNode = get_json_node(node, "parameters");
        json_node* scaleNode = get_json_node(parametersNode, "scale");
        json_node* shapeNode = get_json_node(parametersNode, "shape");

        if(scaleNode == NULL || shapeNode == NULL){
            fprintf(stderr, "Weibull distribution should be parametrized with scale and shape\n");
            exit(13);
        }
        scale = distmodel_parse_parameter(scaleNode, hash, "", 0, INFINITY);
        shape = distmodel_parse_parameter(shapeNode, hash, "", 0, INFINITY);
    }

    Parameters_move(parameters, scale);
    Parameters_move(parameters, shape);

    DistributionModel* dm = new_WeibullDistributionModel_with_parameters(parameters, x);

    dm->shift = get_json_node_value_double(node, "shift", dm->shift);

    Model* model = new_DistributionModel2(id, dm);

    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");

    free_Parameters(x);
    free_Parameters(parameters);

    return model;
}
