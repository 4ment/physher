//
//  ophmc.c
//  physher
//
//  Created by Mathieu Fourment on 11/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "ophmc.h"

#include <math.h>
#include <string.h>

#include "matrix.h"
#include "utilsgsl.h"
#include "gaussian.h"

// op->parameters layout for the HMC operator:
//   [0] step size epsilon (the value used by the next proposal)
//   [1] number of leapfrog steps
//   [2] mu          - anchor point log(10*epsilon0), set on first adaptation
//   [3] Hbar        - running average of (target - accept_prob)
//   [4] logEpsBar   - smoothed (averaged) log step size
//   [5] m           - adaptation step counter (0 => not yet initialized)
#define HMC_PARAM_COUNT 6
#define HMC_MU 2
#define HMC_HBAR 3
#define HMC_LOGEPSBAR 4
#define HMC_M 5

// HMC step-size adaptation via Nesterov dual averaging (Hoffman & Gelman 2014,
// "The No-U-Turn Sampler", Algorithm 5). It drives the mean Metropolis
// acceptance probability towards op->target. Adaptation runs for the whole
// chain: the per-step change vanishes (sqrt(m)/gamma against Hbar -> 0, and the
// m^-kappa smoothing), so the diminishing-adaptation condition holds and the
// stationary distribution is preserved without a hard warmup cut-off.
void operator_hmc_optimize(Operator* op, double logAlpha){
	long count, accepted;
	// Respect the tuning delay; we don't use the cumulative-rate stats here.
	if(!operator_tuning_stats(op, &count, &accepted)) return;

	static const double gamma = 0.05;
	static const double t0 = 10.0;
	static const double kappa = 0.75;

	// First call past the delay: anchor mu to 10x the initial step size, the
	// canonical dual-averaging initialization (start high, shrink towards it).
	if (op->parameters[HMC_M] == 0.0) {
		op->parameters[HMC_MU] = log(10.0 * op->parameters[0]);
		op->parameters[HMC_HBAR] = 0.0;
		op->parameters[HMC_LOGEPSBAR] = 0.0;
	}
	double m = (op->parameters[HMC_M] += 1.0);

	// Current proposal's acceptance probability, clamped to [0,1]. A divergent
	// trajectory gives 0, which pulls the step size down. NaN must be tested
	// before fmin: C99 fmin(0, NaN) returns 0, so exp(fmin(0, NaN)) would be 1
	// (a "perfect acceptance") and drive the step size up — a runaway.
	double alpha = isnan(logAlpha) ? 0.0 : exp(fmin(0.0, logAlpha));

	double etaH = 1.0 / (m + t0);
	op->parameters[HMC_HBAR] =
		(1.0 - etaH) * op->parameters[HMC_HBAR] + etaH * (op->target - alpha);

	double logEps = op->parameters[HMC_MU] - sqrt(m) / gamma * op->parameters[HMC_HBAR];
	double etaX = pow(m, -kappa);
	op->parameters[HMC_LOGEPSBAR] =
		etaX * logEps + (1.0 - etaX) * op->parameters[HMC_LOGEPSBAR];

	// Propose with the raw epsilon during adaptation so Hbar gets matched
	// feedback; logEpsBar is the smoothed value to freeze at if warmup ends.
	op->parameters[0] = exp(logEps);
}

bool operator_hmc(Operator* op, double* logHR){
	size_t paramCount = Parameters_count(op->x);
	size_t dim = Parameters_size(op->x);
	Model* posterior = op->models[0];
	double stepSize = op->parameters[0];
	int steps = op->parameters[1];

	double* momentum0 = dvector(dim);
	double* momentum = dvector(dim);
	double* position = dvector(dim);

	for(size_t i = 0; i < dim; i++){
		momentum0[i] = rnorm();
	}
	memcpy(momentum, momentum0, dim*sizeof(double));
	size_t offset = 0;
	for(size_t i = 0; i < paramCount; i++){
		Parameter* p = Parameters_at(op->x, i);
		const double* values = Parameter_values(p);
		memcpy(position + offset, values, Parameter_size(p)*sizeof(double));
		offset += Parameter_size(p);
	}

	// const double U0 = -posterior->logP(posterior);
	Parameters_zero_grad(op->x);
	posterior->gradient(posterior, op->x);
	offset = 0;
	for(size_t i = 0; i < paramCount; i++){
		Parameter* p = Parameters_at(op->x, i);
		for(size_t j = 0; j < Parameter_size(p); j++){
			double dU = -p->grad[j];
			momentum[offset] -= stepSize/2.0 * dU;
			offset++;
		}
	}

	// A trajectory diverges when a leapfrog step overshoots into a region where
	// the posterior (and hence its gradient) is non-finite. Continuing to
	// integrate would evaluate gradients on a degenerate tree; instead we abort
	// and reject (logHR = -inf), which also lets the step-size tuner shrink it.
	bool diverged = false;
	for (size_t s = 0; s < steps && !diverged; s++) {
		for(size_t i = 0; i < dim; i++){
			position[i] += stepSize * momentum[i];
		}
		Parameters_set_values(op->x, position);

		Parameters_zero_grad(op->x);
		posterior->gradient(posterior, op->x);
		offset = 0;
		for(size_t i = 0; i < paramCount; i++){
			Parameter* p = Parameters_at(op->x, i);
			for(size_t j = 0; j < Parameter_size(p); j++){
				double dU = -p->grad[j];
				if (!isfinite(dU)) {
					diverged = true;
					break;
				}
				momentum[offset] -= stepSize * dU;
				offset++;
			}
			if (diverged) break;
		}
	}

	if (diverged) {
		free(momentum0);
		free(momentum);
		free(position);
		*logHR = -INFINITY;
		return true;
	}

	offset = 0;
	for(size_t i = 0; i < paramCount; i++){
		Parameter* p = Parameters_at(op->x, i);
		for(size_t j = 0; j < Parameter_size(p); j++){
			double dU = -p->grad[j];
			momentum[offset] += stepSize/2.0 * dU;
			offset++;
		}
	}

	double K0 = 0.0;
	double K1 = 0.0;
    for(size_t i = 0; i < dim; i++){
        K0 += momentum0[i] * momentum0[i];
        K1 += momentum[i] * momentum[i];
    }
    K0 /= 2.0;
    K1 /= 2.0;
	// const double U1 = -posterior->logP(posterior);
	// *logHR = (U0 + K0) - (U1 + K1);
	*logHR = K0 - K1;

	free(momentum0);
	free(momentum);
	free(position);

	return true;
}

Operator* new_HMCOperator_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"algorithm",
		"coalescent",
		"delay",
		"model",
		"parameters",
		"stepsize",
		"steps",
		"target",
		"tree",
		"weight",
		"x"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	Operator* op = malloc(sizeof(Operator));
	const char* id_string = get_json_node_value_string(node, "id");
	op->weight = get_json_node_value_double(node, "weight", 1);
	op->name = String_clone(id_string);
	
	op->x = new_Parameters(1);
	get_parameters_references2(node, hash, op->x, "x");
	char* ref = get_json_node_value_string(node, "model");
	op->models = malloc(sizeof(Model*));
	// posterior model
	op->models[0] = Hashtable_get(hash, ref+1);
	op->models[0]->ref_count++;
	op->model_count = 1;
	
	op->propose = operator_hmc;
	op->optimize = operator_hmc_optimize;
	// dvector() zero-fills, so the dual-averaging state ([2..5], in particular
	// the m counter [5]) starts at 0 and is initialized on the first tune call.
	op->parameters = dvector(HMC_PARAM_COUNT);
	op->parameters[0] = get_json_node_value_double(node, "stepsize", 0.1);
	op->parameters[1] = get_json_node_value_double(node, "steps", 5);
	
	op->rejected_count = 0;
	op->accepted_count = 0;
	op->failure_count = 0;
	op->tuning_delay = get_json_node_value_size_t(node, "delay", 0);
	op->accepted_at_delay = 0;
	op->count_at_delay = 0;
	op->tuning_started = false;
	op->all = false;
	op->target = get_json_node_value_double(node, "target", 0.8);
	op->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	return op;
}
