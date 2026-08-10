// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "optimizer.h"

#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <strings.h>

#include "parameters.h"
#include "utils.h"
#include "brent.h"
#include "em.h"
#include "matrix.h"
#include "powell.h"
#include "bfgs.h"
#include "frpmrn.h"
#include "gradascent.h"
#include "topologyopt.h"
#include "tree.h"
#include "treelikelihood.h"
#include "tracelogger.h"
#include "modelfactory.h"

static double model_logP( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;
	double logP = model->logP(model);
	if (grad != NULL) {
		for(size_t i = 0; i < Parameters_count(params); i++){
			Parameter_zero_grad(Parameters_at(params, i));
		}
		model->gradient(model, params);
		
		size_t index = 0;
		for(size_t i = 0; i < Parameters_count(params); i++){
			Parameter* p = Parameters_at(params, i);
			for(size_t j = 0; j < Parameter_size(p); j++){
				grad[index] = p->grad[j];
				index++;
			}
		}
	}
	//	printf("%f\n", logP);
	return logP;
}

double model_negative_logP( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;
	double logP = model->logP(model);
//		printf("%f\n", model->logP(model));
	if (grad != NULL) {
		for(size_t i = 0; i < Parameters_count(params); i++){
			Parameter_zero_grad(Parameters_at(params, i));
		}
		model->gradient(model, params);
		
		size_t index = 0;
		for(size_t i = 0; i < Parameters_count(params); i++){
			Parameter* p = Parameters_at(params, i);
			for(size_t j = 0; j < Parameter_size(p); j++){
				grad[index] = p->grad[j];
				index++;
			}
		}
	}
	return -logP;
}

static void _gradient( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;
	for(size_t i = 0; i < Parameters_count(params); i++){
		Parameter_zero_grad(Parameters_at(params, i));
	}
	model->gradient(model, params);

	if(grad != NULL){
		size_t index = 0;
		for(size_t i = 0; i < Parameters_count(params); i++){
			Parameter* p = Parameters_at(params, i);
			for(size_t j = 0; j < Parameter_size(p); j++){
				grad[index] = p->grad[j];
				index++;
			}
		}
	}
}

static void _negative_gradient( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;
	for(size_t i = 0; i < Parameters_count(params); i++){
		Parameter_zero_grad(Parameters_at(params, i));
	}
	model->gradient(model, params);
	size_t index = 0;
	for(size_t i = 0; i < Parameters_count(params); i++){
		Parameter* p = Parameters_at(params, i);
		for(size_t j = 0; j < Parameter_size(p); j++){
			grad[index] = -p->grad[j];
			index++;
		}
	}
}

static void _reset(void* data){
	Model* model = (Model*)data;
	// model->reset(model);
}

static bool dummy_update_data( void *data, Parameters *p){return false;}

static bool xStop( const Parameters *x,  double *xold, const double tolx);

struct _Optimizer{
	opt_algorithm algorithm;
	unsigned int dimension;
	
	opt_func f;
    opt_grad_func grad_f;
	void (*reset)(void*);
	void *data;
//    Model* model;
	Parameters* parameters;
	Model* treelikelihood;

	OptStopCriterion stop;
	
	opt_update_data update;
	int verbosity;
	OptimizerSchedule* schedule;
    
	char* checkpoint_file;
	int checkpoint_frequency;
    // for stochastic gradient
    double* etas;
	size_t eta_count;
    bool ascent;
	size_t threads;
	bool maximize;
	Trace* logger;
};

opt_result topology_optimize(TopologyOptimizer* topopt, double *fmin){
	*fmin = -topopt->optimize(topopt);
	return OPT_SUCCESS;
}

// EM has no parameter list and no objective of its own: it reads the free-rate
// site model off `opt->treelikelihood` and updates it in place (see em.h). The
// value reported back is the *target*'s, evaluated through opt->f so the sign
// convention matches every other entry -- which also means a target carrying
// priors is scored, but not maximized, by this step: EM maximizes the likelihood
// alone.
static opt_result em_optimize( Optimizer *opt, double *fmin ){
	if (opt->treelikelihood == NULL) {
		fprintf(stderr, "EM has no tree likelihood to read the mixture off\n");
		return OPT_ERROR;
	}
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)opt->treelikelihood->obj;
	SiteModelEM em = SiteModel_optimize_freerate_EM(tlk, opt->stop.iter_max,
	                                                opt->stop.tolx);
	*fmin = opt->f(NULL, NULL, opt->data);
	opt->stop.iter = em.steps;
	// What the per-class M-steps spent, plus one E-step traversal per step. Both
	// are full traversals, so the sum is comparable with what the other entries
	// of a meta schedule report.
	opt->stop.f_eval_current += em.evaluations + em.steps;
	if (isnan(em.logP)) return OPT_ERROR;
	return em.converged ? OPT_SUCCESS : OPT_MAXITER;
}

opt_result serial_brent_optimize_tree( Model* mtlk, opt_func f, void *data, OptStopCriterion *stop, double *fmin ){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)mtlk->obj;
	Tree* tree = tlk->tree;
	Node** nodes = Tree_get_nodes(tree, POSTORDER);
	tlk->node_upper = NULL;
	tlk->use_upper = true;
	tlk->update_upper = true;

	// initialize lower and upper
	SingleTreeLikelihood_update_uppers(tlk);
	//for(int j = 0; j < stop->iter_min; j++)
	for(int i = 0; i < Tree_node_count(tree); i++){
		Node* node = nodes[i];
		// skips the root and the child of the root pinned to a zero branch
		if(!Node_has_distance(node)) continue;
		if(tlk->node_upper == NULL) tlk->node_upper = node;

		// iter is Brent's own inner-loop counter and restarts per branch. The
		// evaluation count and the clock must not: they are budgets over the
		// whole call, and restarting them here gave every branch a fresh
		// allowance and reported only the last branch's work.
		stop->iter = 0;
		stop->count = 0;
#ifdef UPPER_PARTIALS
		printf("brent\n");
#endif
		// printf("%s %f %f %f\n", node->distance->name, node->distance->value[node->id], Constraint_lower(node->distance->cnstr), Constraint_upper(node->distance->cnstr));
		opt_result status = brent_optimize2(node->distance, node->branch_index, f, data, stop, fmin);
		// printf(" %d %f %f\n",  node->id, node->distance->value[node->id], *fmin);
#ifdef UPPER_PARTIALS
		printf("%f %s\n", -*fmin, nodes[i]->name);
#endif
		//tlk->node_upper = nodes[i];
	}
    
	tlk->use_upper = false;
	return OPT_SUCCESS;
}

// Judge one step of an optimizer: a sweep of a meta schedule, a Powell direction
// set, a conjugate-gradient iteration. `before` and `after` are the objective at
// the start and end of the step; it is a negative log-likelihood, so smaller is
// better and `before - after` is the gain. The caller decides what `before` is --
// meta passes the value at the head of the sweep, opt_check_stop passes the value
// at the last check that made progress.
//
// Three things distinguish this from the `before - after < tolfx` test it
// replaces:
//
//  - `iter_min` (the JSON `min` key) is honoured. It was parsed and then never
//    read by anything, so a single flat sweep ended the run. That is exactly the
//    case a schedule needs to push through: sweeps can crawl along a curved
//    ridge before an optimizer further down the list breaks out of it.
//  - `patience` requires several consecutive flat steps rather than one.
//  - a step that ends *higher* than it started is reported as OPT_PROGRESS_WORSE.
//    The
//    old test folded it into the "gain below tolerance" branch, because a
//    negative gain is also less than tolfx, and returned OPT_SUCCESS from a
//    point worse than the one the sweep started at.
//
// The tolerance stays absolute. A log-likelihood difference carries the same
// meaning whatever the total, so `tolfx` is a number of log units, not a
// fraction of the objective -- scaling it by the magnitude would ask a run whose
// likelihood is -1e6 to accept a gain of tens of units as convergence. The only
// scaling applied is a floor at the representable resolution of `after`, so a
// tolerance smaller than an ulp cannot demand a difference doubles cannot express.
//
// `stop->stall` carries the number of consecutive steps whose gain was below
// tolerance; the caller zeroes it before the first step.
opt_progress opt_check_progress( OptStopCriterion *stop, double before, double after ){
	double tol = fmax(stop->tolfx, 8.0*DBL_EPSILON*fabs(after));
	double gain = before - after;

	if (gain < -tol) {
		stop->stall = 0;
		return OPT_PROGRESS_WORSE;
	}
	if (gain > tol) {
		stop->stall = 0;
		return OPT_PROGRESS_ONGOING;
	}
	stop->stall++;
	if (stop->stall < stop->patience || stop->iter < stop->iter_min) {
		return OPT_PROGRESS_ONGOING;
	}
	return OPT_PROGRESS_CONVERGED;
}

// Sweeps over which a schedule entry's configured `rounds` are honoured. The
// extra rounds only pay off while the parameters are far from the optimum;
// after this many sweeps every entry drops to one round.
#define META_MULTIROUND_SWEEPS 2


// ---------------------------------------------------------------------------
// Progress table for the meta optimizer
//
// One row per schedule entry per sweep, then a total row closing the sweep.
// The point of the table is to show where the gain in the objective comes from
// and what each entry costs, so every row carries the value the entry left the
// objective at, what it gained, and how much work that took. Columns are
// plain ASCII and fixed width, like the operator summary in mcmc.c.
// ---------------------------------------------------------------------------

#define META_NAME_WIDTH_MAX 24
#define META_ALGO_WIDTH 11  // "BRENTSERIAL"

// Label a schedule entry. The topology test comes first: a topology entry can
// also carry a tree likelihood, and it does not optimize branch lengths.
static const char* meta_entry_name(const Optimizer* opt){
	if (opt->algorithm == OPT_TOPOLOGY) return "topology";
	if (opt->algorithm == OPT_EM) return "free rates";
	// Every other algorithm holding a tree likelihood is there to sweep its
	// branch lengths.
	if (opt->treelikelihood != NULL) return "branches";
	if (opt->parameters != NULL && Parameters_count(opt->parameters) > 0){
		const char* group = Parameters_name2(opt->parameters);
		if (group != NULL) return group;
		return Parameters_name(opt->parameters, 0);
	}
	return "(unnamed)";
}

// Whatever an entry can say about itself beyond the objective: the number of
// moves the topology optimizer accepted, the value the parameters landed on.
// `name` is the label already in the Component column, so a single parameter
// named after its own entry contributes its value alone.
static void meta_entry_detail(const Optimizer* opt, const char* name, char* buffer, size_t size){
	if (opt->algorithm == OPT_TOPOLOGY){
		TopologyOptimizer* topopt = opt->data;
		snprintf(buffer, size, "%d move%s", topopt->moves, topopt->moves == 1 ? "" : "s");
	}
	else if (opt->parameters != NULL && Parameters_count(opt->parameters) > 0){
		size_t count = Parameters_count(opt->parameters);
		const char* first = Parameters_name(opt->parameters, 0);
		double value = Parameters_value(opt->parameters, 0);
		if (count == 1 && strcmp(first, name) == 0){
			snprintf(buffer, size, "%.6g", value);
		}
		else if (count == 1){
			snprintf(buffer, size, "%s = %.6g", first, value);
		}
		else{
			snprintf(buffer, size, "%zu parameters, %s = %.6g", count, first, value);
		}
	}
	else{
		buffer[0] = '\0';
	}
}

static int meta_name_width(const OptimizerSchedule* schedule){
	size_t width = strlen("Component");
	for (int i = 0; i < schedule->count; i++){
		size_t len = strlen(meta_entry_name(schedule->optimizers[i]));
		if (len > width) width = len;
	}
	if (width > META_NAME_WIDTH_MAX) width = META_NAME_WIDTH_MAX;
	return (int)width;
}

// Right-aligned in the Elapsed column, in whichever unit keeps the number
// readable: 0.42s, 12.3s, 4m18s.
static void meta_format_duration(double seconds, char* buffer, size_t size){
	if (seconds < 60.0){
		snprintf(buffer, size, "%.2fs", seconds);
	}
	else if (seconds < 3600.0){
		snprintf(buffer, size, "%dm%02ds", (int)(seconds/60), (int)fmod(seconds, 60.0));
	}
	else{
		snprintf(buffer, size, "%dh%02dm", (int)(seconds/3600), (int)fmod(seconds/60, 60.0));
	}
}

// How the run ended, for the closing line of the table.
static const char* meta_outcome(opt_result result){
	switch (result) {
		case OPT_SUCCESS: return "Converged";
		case OPT_MAXITER: return "Stopped on the sweep limit";
		case OPT_MAXEVAL: return "Stopped on the evaluation limit";
		case OPT_MAXTIME: return "Stopped on the time limit";
		case OPT_FAIL:    return "Ended worse than it started";
		default:          return "Ended";
	}
}

static void meta_print_row(int name_width, size_t sweep, const char* name, const char* algorithm,
                           double logL, double gain, size_t evals, double seconds, const char* detail){
	char duration[32];
	meta_format_duration(seconds, duration, sizeof(duration));
	printf("  %5zu  %-*.*s  %-*s  %15.4f  %+12.4f  %8zu  %8s", sweep, name_width, name_width,
	       name, META_ALGO_WIDTH, algorithm, logL, gain, evals, duration);
	if (detail != NULL && detail[0] != '\0') printf("  %s", detail);
	putchar('\n');
}

// Spans every column except the free-width Detail one: the fixed columns are
// 5 + name_width + 11 + 15 + 12 + 8 + 8 wide with a two-space gutter between them.
static void meta_print_rule(int name_width){
	printf("  ");
	for (int i = 0; i < 71 + name_width; i++) putchar('-');
	putchar('\n');
}

// `logL` is the objective before any entry of the schedule has run: sweep 0 of
// the table, so every Delta below it can be read against the point the run
// started from. It costs the one evaluation already charged to the run.
static void meta_print_header(const OptimizerSchedule* schedule, const OptStopCriterion* stop,
                              int name_width, bool maximize, double logL){
	printf("\nMeta optimizer: %d component%s, up to %zu sweep%s, tolerance %g\n",
	       schedule->count, schedule->count == 1 ? "" : "s",
	       stop->iter_max, stop->iter_max == 1 ? "" : "s", stop->tolfx);
	printf("  %5s  %-*s  %-*s  %15s  %12s  %8s  %8s  %s\n", "Sweep", name_width, "Component",
	       META_ALGO_WIDTH, "Algorithm", maximize ? "logL" : "f", "Delta", "Evals", "Elapsed",
	       "Detail");
	meta_print_rule(name_width);
	printf("  %5zu  %-*s  %-*s  %15.4f  %12s  %8zu\n", (size_t)0, name_width, "start",
	       META_ALGO_WIDTH, "", logL, "-", (size_t)1);
}

static opt_result meta_optimize( Optimizer* opt_meta, double *fmin ){
	opt_func f = opt_meta->f;
	void* data = opt_meta->data;
	OptStopCriterion* stop = &opt_meta->stop;
	OptimizerSchedule* schedule = opt_meta->schedule;
	const int verbosity = opt_meta->verbosity;
	// The objective is the negative log-likelihood when the run maximizes, so
	// every value the table reports is negated back to a log-likelihood. A
	// minimizing run reports the objective itself.
	const double sign = opt_meta->maximize ? -1.0 : 1.0;
	const int name_width = meta_name_width(schedule);
	struct timespec time_start;
	time_monotonic(&time_start);

	double lnl = f(NULL, NULL, data);
	const double lnl_start = lnl;
	double fret = lnl;
	*fmin = lnl;
	stop->stall = 0;
	stop->iter = 0;
	stop->f_eval_current = 1;
	opt_result result = OPT_MAXITER;

	if (verbosity > 0) meta_print_header(schedule, stop, name_width, opt_meta->maximize, sign*lnl_start);

	for (size_t sweep = 0; sweep < stop->iter_max; sweep++) {
		double lnl_current = lnl;
		const size_t evals_before = stop->f_eval_current;
		struct timespec sweep_start;
		time_monotonic(&sweep_start);
		for (int i = 0; i < schedule->count; i++) {
			Optimizer* opt = schedule->optimizers[i];
			const double lnl_before = lnl;
			size_t entry_evals = 0;
			bool entry_failed = false;
			struct timespec entry_start;
			time_monotonic(&entry_start);
			// Read the round count, do not overwrite it. This used to assign
			// schedule->rounds[i] = 1 on the third sweep, which permanently
			// rewrote the schedule: an Optimizer reused across bootstrap
			// replicates ran the configured rounds for replicate 0 and a single
			// round for every replicate after it.
			int rounds = (sweep < META_MULTIROUND_SWEEPS) ? schedule->rounds[i] : 1;
			double local_fret;
			for (int k = 0; k < rounds; k++){
				local_fret = fret;
				opt_result status;
				// Entries report how many times they evaluated the objective;
				// zero the counter first so what comes back is this call's work
				// and not an accumulation over every sweep so far.
				// serial_brent_optimize_tree does not go through opt_optimize, so
				// it would otherwise never be reset.
				opt->stop.f_eval_current = 0;
				// A tree likelihood on an entry means "sweep its branch lengths"
				// -- except for the two algorithms that hold one for their own
				// reasons and dispatch through opt_optimize like everything else.
				if(opt->treelikelihood != NULL && opt->algorithm != OPT_TOPOLOGY &&
				   opt->algorithm != OPT_EM){
					status = serial_brent_optimize_tree(opt->treelikelihood, opt->f, opt->data, &opt->stop, &fret);
				}
				else{
					status = opt_optimize( opt, &fret);
				}
				stop->f_eval_current += opt->stop.f_eval_current;
				entry_evals += opt->stop.f_eval_current;
				// A failed child has not necessarily written fret -- scaler_optimize
				// returns OPT_FAIL without touching it -- so the sweep would carry
				// the previous entry's value forward as if this one had run. Resync
				// on the target instead of trusting the report. The run itself is
				// not abandoned: Brent legitimately reports OPT_FAIL when it rolls
				// back to its incoming value.
				if (status == OPT_FAIL || status == OPT_ERROR) {
					fret = f(NULL, NULL, data);
					stop->f_eval_current++;
					entry_evals++;
					entry_failed = true;
				}
			}
			lnl = fret;

			if (verbosity > 0) {
				char detail[128];
				const char* name = meta_entry_name(opt);
				meta_entry_detail(opt, name, detail, sizeof(detail));
				if (entry_failed) {
					// Reported in the row rather than on stderr so it stays next
					// to the sweep it happened in.
					size_t used = strlen(detail);
					snprintf(detail + used, sizeof(detail) - used, "%s[failed]",
					         used > 0 ? "  " : "");
				}
				meta_print_row(name_width, sweep + 1, name,
				               OPT_ALGORITHMS[opt->algorithm], sign*lnl, sign*(lnl - lnl_before),
				               entry_evals, time_elapsed(&entry_start), detail);
			}
		}

		// Rescaling is switched on inside the likelihood when the partials
		// underflow; turn it back off once a sweep so the cheaper unscaled path
		// is used again whenever the parameters allow it. Any entry of the
		// schedule may be the one holding the tree -- this used to look only at
		// optimizers[0], so a schedule that did not list the branch-length
		// optimizer first never reached it.
		for (int i = 0; i < schedule->count; i++) {
			if (schedule->optimizers[i]->treelikelihood == NULL) continue;
			SingleTreeLikelihood* tlk = schedule->optimizers[i]->treelikelihood->obj;
			if (tlk->scale) {
				SingleTreeLikelihood_use_rescaling(tlk, false);
				SingleTreeLikelihood_update_all_nodes( tlk );
			}
		}

		// Re-evaluate the target rather than trusting whatever the last entry of
		// the schedule reported. Serial Brent reports an upper-partial likelihood,
		// the topology optimizer reports its own accounting, and the rescaling
		// toggle above has just invalidated every partial -- so the two numbers
		// the criterion compares could otherwise come from different code paths.
		// One full evaluation per sweep is negligible next to the hundreds the
		// schedule itself performs.
		lnl = fret = f(NULL, NULL, data);
		stop->f_eval_current++;
		stop->iter = sweep + 1;
		*fmin = lnl;

		opt_progress convergence = opt_check_progress(stop, lnl_current, lnl);

		if (verbosity > 0) {
			char detail[64] = "";
			if (convergence == OPT_PROGRESS_CONVERGED) {
				snprintf(detail, sizeof(detail), "converged");
			}
			else if (convergence == OPT_PROGRESS_WORSE) {
				snprintf(detail, sizeof(detail), "worse than the previous sweep");
			}
			else if (stop->stall > 0) {
				// Flat sweeps so far against the number `patience` demands before
				// the run is called converged.
				snprintf(detail, sizeof(detail), "flat %zu/%zu", stop->stall, stop->patience);
			}
			meta_print_row(name_width, sweep + 1, "= sweep", "", sign*lnl,
			               sign*(lnl - lnl_current), stop->f_eval_current - evals_before,
			               time_elapsed(&sweep_start), detail);
			putchar('\n');
		}

		if (convergence == OPT_PROGRESS_CONVERGED) {
			result = OPT_SUCCESS;
			break;
		}

		// Budgets are enforced at sweep granularity: an entry of the schedule
		// cannot be interrupted once it is running, so a limit is honoured to
		// within one sweep. Checked after the convergence test so a run that
		// finishes on the same sweep it runs out of budget still reports success.
		opt_result limit = opt_check_limits(stop);
		if (limit != OPT_KEEP_GOING) {
			result = limit;
			break;
		}
	}

	// Converging is not the same as improving. The schedule has no way to roll
	// the model back -- its entries own disjoint parameter sets and the
	// branch-length entry owns none -- so the best it can do is refuse to call a
	// net regression a success.
	if (result == OPT_SUCCESS && lnl > lnl_start + stop->tolfx * fmax(1.0, fabs(lnl_start))) {
		result = OPT_FAIL;
	}

	if (verbosity > 0) {
		char duration[32];
		meta_format_duration(time_elapsed(&time_start), duration, sizeof(duration));
		meta_print_rule(name_width);
		printf("  %s after %zu sweep%s: %s %.4f -> %.4f (%+.4f) in %zu evaluations, %s\n\n",
		       meta_outcome(result), stop->iter, stop->iter == 1 ? "" : "s",
		       opt_meta->maximize ? "logL" : "f", sign*lnl_start, sign*lnl,
		       sign*(lnl - lnl_start), stop->f_eval_current, duration);
	}
	return result;
}



Optimizer * new_Optimizer( opt_algorithm algorithm ) {
	Optimizer * opt;
	opt = (Optimizer*) malloc(sizeof(struct _Optimizer));
	assert(opt);
	opt->algorithm = algorithm;
	opt->f = NULL;
    opt->grad_f = NULL;
	opt->data = NULL;
	opt->parameters = NULL;
	opt->treelikelihood = NULL;
	opt->dimension = 0;
	
	opt->stop.iter_min = 1;
	opt->stop.iter_max = 200;
	opt->stop.iter = 0;

	opt->stop.time_max = 0;
	opt->stop.time_start = 0;
	opt->stop.time_end = 0;
	opt->stop.time_current = 0;
	
	opt->stop.f_eval_max = 0;
	opt->stop.f_eval_current = 0;
	
	opt->stop.tolfx = 0;
	opt->stop.tolx = OPT_XTOL;
	opt->stop.tolg = 1.e-5;
	opt->stop.patience = 1;
	opt->stop.stall = 0;

	opt->stop.oldx = NULL;
	opt->stop.oldfx = 0;
	
	opt->stop.count = 0;
	opt->verbosity = 1;
	opt->update = dummy_update_data;
	opt->schedule = NULL;
    
	opt->stop.frequency_check = 100;
    opt->etas = NULL;
	opt->eta_count = 0;
    opt->ascent = true;
	opt->reset = NULL;

	opt->checkpoint_file = NULL;
	opt->checkpoint_frequency = 0;
	opt->logger = NULL;
	return opt;
}

// need to change:
// data treelikelihood
// f: maybe different log likleihood for tree with upper
Optimizer* clone_Optimizer(Optimizer *opt, void* data, Parameters* parameters){
	Optimizer* clone = (Optimizer*) malloc(sizeof(struct _Optimizer));
	assert(clone);
	clone->algorithm = opt->algorithm;
	clone->f = opt->f;
	clone->data = data;
	clone->dimension = opt->dimension;
	
	clone->stop.iter_min = opt->stop.iter_min;
	clone->stop.iter_max = opt->stop.iter_max;
	clone->stop.iter = opt->stop.iter;
	clone->stop.frequency_check = opt->stop.frequency_check;
	
	clone->stop.time_max = opt->stop.time_max;
	clone->stop.time_start = opt->stop.time_start;
	clone->stop.time_end = opt->stop.time_end;
	clone->stop.time_current = opt->stop.time_current;
	
	clone->stop.f_eval_max = opt->stop.f_eval_max;
	clone->stop.f_eval_current = opt->stop.f_eval_current;
	
	clone->stop.tolfx = opt->stop.tolfx;
	clone->stop.tolx = opt->stop.tolx;
	clone->stop.tolg = opt->stop.tolg;
	clone->stop.patience = opt->stop.patience;
	clone->stop.stall = opt->stop.stall;

	clone->stop.oldx = opt->stop.oldx;
	clone->stop.oldfx = opt->stop.oldfx;
	
	clone->stop.count = opt->stop.count;
	clone->verbosity = opt->verbosity;
	clone->update = opt->update;
	clone->schedule = NULL;
	clone->parameters = NULL;
	clone->treelikelihood = NULL;
    
    clone->ascent = opt->ascent;
	clone->etas = clone_dvector(opt->etas, opt->eta_count);
	clone->eta_count = opt->eta_count;
    clone->grad_f = opt->grad_f;
	
	if(opt->schedule != NULL){
		clone->schedule = (OptimizerSchedule*)malloc(sizeof(OptimizerSchedule));
		clone->schedule->capacity = opt->schedule->capacity;
		clone->schedule->count = opt->schedule->count;
		clone->schedule->optimizers = (Optimizer**)calloc(opt->schedule->capacity, sizeof(Optimizer*));
		clone->schedule->rounds = ivector(opt->schedule->capacity);
		memcpy(clone->schedule->rounds, opt->schedule->rounds, sizeof(int)*opt->schedule->count);
		
		for (int i = 0; i < clone->schedule->count; i++) {
			clone->schedule->optimizers[i] = clone_Optimizer(opt->schedule->optimizers[i], data, parameters);
			
			if(Parameters_count(opt->parameters) > 0){
				opt->parameters = new_Parameters(Parameters_count(opt->parameters));
			}
			
			// Find matching parameters
			for (int j = 0; j < Parameters_count(opt->parameters); j++) {
				int k = 0;
				for (k = 0; k < Parameters_count(parameters); k++) {
					if(strcmp(Parameters_name(opt->parameters, j), Parameters_name(parameters, k)) == 0){
						Parameters_add(clone->parameters, Parameters_at(parameters, k));
						break;
					}
				}
				if(k == Parameters_count(parameters)){
					printf("not found %s\n", Parameters_name(parameters, k));
					for (k = 0; k < Parameters_count(parameters); k++) {
						printf("%s\n", Parameters_name(parameters, k));
					}
					exit(1);
				}
			}
		}
	}
	clone->threads = opt->threads;
	return clone;
}

void free_Optimizer( Optimizer *opt ){
	if(opt->algorithm == OPT_TOPOLOGY){
		TopologyOptimizer* topopt = opt->data;
		free_TopologyOptimizer(topopt);
	}
	opt->data = NULL;
	free_Parameters(opt->parameters);
	if(opt->etas != NULL){
		free(opt->etas);
	}
	//if(opt->treelikelihood != NULL) opt->treelikelihood->free(opt->treelikelihood);
	if(opt->schedule != NULL){
		for (int i = 0; i < opt->schedule->count; i++) {
			free_Optimizer(opt->schedule->optimizers[i]);
		}
		free(opt->schedule->optimizers);
		free(opt->schedule->rounds);
		free(opt->schedule);
	}
	if(opt->checkpoint_file != NULL){
		free(opt->checkpoint_file);
	}
	if(opt->logger != NULL){
		opt->logger->free(opt->logger);
	}
	free(opt);
}

void opt_set_objective_function( Optimizer *opt, opt_func f ){
	if( opt != NULL ){
		opt->f = f;
	}
}

void opt_set_data( Optimizer *opt, void *data ){
	if( opt != NULL ){
		opt->data = data;
	}
}

void opt_set_parameters( Optimizer *opt, const Parameters *parameters ){
	if( opt != NULL ){
		if(opt->parameters == NULL){
			opt->parameters = new_Parameters(Parameters_count(parameters));
			if(Parameters_name2(parameters) != NULL){
				Parameters_set_name2(opt->parameters, Parameters_name2(parameters));
			}
		}
		else{
			Parameters_removeAll(opt->parameters);
		}
		Parameters_add_parameters(opt->parameters, parameters);
	}
}

void opt_set_treelikelihood( Optimizer *opt, Model* treelikelihood){
	opt->treelikelihood = treelikelihood;
}

void opt_set_max_evaluation( Optimizer *opt, const size_t maxeval ){
	if( opt != NULL ){
		opt->stop.f_eval_max = maxeval;
	}
}

void opt_set_max_iteration( Optimizer *opt, const size_t maxiter ){
	if( opt != NULL ){
		opt->stop.iter_max = maxiter;
	}
}

void opt_set_min_iteration( Optimizer *opt, const size_t miniter ){
	if( opt != NULL ){
		opt->stop.iter_min = miniter;
	}
}

void opt_set_patience( Optimizer *opt, const size_t patience ){
	if( opt != NULL ){
		opt->stop.patience = patience;
	}
}

void opt_set_verbosity( Optimizer *opt, const int verbosity ){
	if( opt != NULL ){
		opt->verbosity = verbosity;
	}
}

// In seconds by default
void opt_set_time_max( Optimizer *opt, const double maxtime ){
	if( opt != NULL ){
		opt->stop.time_max = maxtime;
	}
}

void opt_set_time_max_minutes( Optimizer *opt, const double maxtime ){
	if( opt != NULL ){
		opt->stop.time_max = maxtime*60;
	}
}

void opt_set_time_max_hours( Optimizer *opt, const double maxtime ){
	if( opt != NULL ){
		opt->stop.time_max = maxtime*3600;
	}
}

void opt_set_time_max_relative( Optimizer *opt, const time_t time, const double factor ){
	if( opt != NULL ){
		opt->stop.time_max = time*factor;
	}
}

void opt_set_tolfx( Optimizer *opt, const double tolfx ){
	if( opt != NULL ){
		opt->stop.tolfx = tolfx;
	}
}

void opt_set_tolx( Optimizer *opt, const double tolx ){
	if( opt != NULL ){
		opt->stop.tolx = tolx;
	}
}

double opt_tolx( Optimizer *opt ){
    return opt->stop.tolx;
}
void opt_set_tolg( Optimizer *opt, const double tolg ){
	if( opt != NULL ){
		opt->stop.tolg = tolg;
	}
}

size_t opt_frequency_check( Optimizer *opt ){
	return opt->stop.frequency_check;
}
void opt_set_frequency_check( Optimizer *opt, const size_t frequency_check ){
	if( opt != NULL ){
		opt->stop.frequency_check = frequency_check;
	}
}

int opt_iterations( Optimizer *opt ){
    return opt->stop.iter;
}

int opt_f_evaluations( Optimizer *opt ){
    return opt->stop.f_eval_current;
}

Parameters* opt_parameters( Optimizer *opt ){
	return opt->parameters;
}
// Meta optimizer

void opt_add_optimizer(Optimizer *opt_meta, Optimizer *opt){
	if(opt_meta->schedule == NULL){
		opt_get_schedule(opt_meta);
	}
	else if(opt_meta->schedule->capacity == opt_meta->schedule->count){
		opt_meta->schedule->capacity++;
		opt_meta->schedule->optimizers = (Optimizer**)realloc(opt_meta->schedule->optimizers, sizeof(Optimizer*)*opt_meta->schedule->capacity);
		opt_meta->schedule->rounds = (int*)realloc(opt_meta->schedule->rounds, sizeof(int)*opt_meta->schedule->capacity);
	}
	opt_meta->schedule->optimizers[opt_meta->schedule->count] = opt;
	opt_meta->schedule->rounds[opt_meta->schedule->count] = 1;
	opt_meta->schedule->count++;
}

OptimizerSchedule* opt_get_schedule(Optimizer *opt_meta){
	if(opt_meta->schedule == NULL){
		opt_meta->schedule = (OptimizerSchedule*)malloc(sizeof(OptimizerSchedule));
		opt_meta->schedule->optimizers = (Optimizer**)malloc(sizeof(Optimizer*));
		opt_meta->schedule->rounds = ivector(1);
		opt_meta->schedule->capacity = 1;
		opt_meta->schedule->count = 0;
	}
	return opt_meta->schedule;
}


// The resource half of opt_check_stop: the wall-clock, iteration and evaluation
// budgets, with no reference to the objective value. Split out so a caller that
// brings its own convergence test -- the meta optimizer -- can honour the
// budgets without also inheriting the value-based test. A zero budget means
// "no limit".
//
// The three are reported in order of urgency rather than being allowed to
// overwrite one another: the old inline version evaluated all three in sequence,
// so an evaluation budget that was still fine erased an already-detected
// timeout and the run kept going.
opt_result opt_check_limits( OptStopCriterion *stop ){
	if( stop->time_max != 0 ){
		time( &stop->time_current );
		if( difftime( stop->time_current, stop->time_start ) > stop->time_max ){
			return OPT_MAXTIME;
		}
	}

	if ( stop->f_eval_max != 0 && stop->f_eval_current + 1 > stop->f_eval_max ) {
		return OPT_MAXEVAL;
	}

	if ( stop->iter_max != 0 && stop->iter > stop->iter_max ) {
		return OPT_MAXITER;
	}

	return OPT_KEEP_GOING;
}

// At least one of these conditions is sufficient to stop the optimization
opt_result opt_check_stop( OptStopCriterion *stop, Parameters *x, double fx ){
	if ( stop->count == 0 ) {
		if(stop->oldx == NULL){
			stop->oldx = dvector(Parameters_size(x));
		}
		Parameters_store_value(x, stop->oldx);
		stop->oldfx = fx;
		stop->count++;
		return OPT_KEEP_GOING;
	}

	opt_result stopflag = opt_check_limits(stop);

//	if( xStop( x, stop->oldx, stop->tolx) ){
//		fprintf(stderr, "tolx criterion: %f\n",stop->tolx);
//		return OPT_SUCCESS;
//	}

	opt_progress progress = opt_check_progress(stop, stop->oldfx, fx);

	// Advance the reference only when the step actually moved the objective --
	// opt_check_progress zeroes the stall counter exactly then. The test this
	// replaces overwrote the reference on every call, which made it "no change
	// since the previous check" rather than "no cumulative progress": a run creeping
	// downhill at just under tolfx per step read as flat at every check while
	// still descending, and raising `patience` did not help because each of the
	// N checks was measured from a fresh reference. Holding it lets the drift
	// accumulate, so a genuine creep is eventually seen as the progress it is.
	if (stop->stall == 0) {
		stop->oldfx = fx;
	}

	if (progress == OPT_PROGRESS_CONVERGED) {
		return OPT_SUCCESS;
	}

	return stopflag;
}

opt_result opt_optimize( Optimizer *opt, double *fmin ){
	if ( opt->stop.time_max != 0 ) {
		time( &opt->stop.time_start );
	}
	opt_result result = OPT_SUCCESS;

	// The parameters to optimize always come from opt->parameters (set via
	// opt_set_parameters); there is no longer a separate per-call parameter list.
	Parameters *ps = opt->parameters;
	opt->dimension = (ps != NULL) ? Parameters_size(ps) : 0;

	// should probably moved somewhere else
	if(opt->dimension > 0){
		opt->stop.oldx = dvector(opt->dimension);
	}
    opt->stop.iter = 0;
    opt->stop.f_eval_current = 0;
    opt->stop.count = 0;
    opt->stop.stall = 0;
	OptimizerCheckpoint checkpointer = {opt->checkpoint_file, opt->checkpoint_frequency};
	
	switch (opt->algorithm ) {
		case OPT_META:{
			result = meta_optimize( opt, fmin );
			break;
		}
		case OPT_POWELL:{
			result = powell_optimize( ps, opt->f, opt->data, &opt->stop, fmin, opt->update );
			break;
		}
		case OPT_BRENT:{
			result = brent_optimize( ps, opt->f, opt->data, &opt->stop, fmin );
			break;
		}
		case OPT_SERIAL_BRENT:{
			if(opt->treelikelihood == NULL){
				result = serial_brent_optimize( ps, opt->f, opt->data, &opt->stop, fmin );
			}
			else{
				result = serial_brent_optimize_tree(opt->treelikelihood, opt->f, opt->data, &opt->stop, fmin);
			}
			break;
		}
		case OPT_BFGS:{
			result = dfpmin_optimize( opt->parameters, opt->f, opt->grad_f, opt->data, &opt->stop, fmin, opt->etas[0]);
			break;
		}
		case OPT_CG_FR:{
			result = frprmn_optimize( opt->parameters, opt->f, opt->grad_f, opt->data, &opt->stop, fmin, OPT_CG_FR );
			break;
        }
        case OPT_CG_PR:{
            result = frprmn_optimize( opt->parameters, opt->f, opt->grad_f, opt->data, &opt->stop, fmin, OPT_CG_PR );
            break;
        }
		case OPT_SG: case OPT_SG_ADAM:{
			double eta = *opt->etas;
			opt_result eta_adapt_result = OPT_SUCCESS;
			if(opt->eta_count > 1){
				eta_adapt_result = optimize_stochastic_gradient_adapt(opt->parameters, opt->f, opt->grad_f, opt->reset, opt->etas, opt->eta_count, opt->data, &opt->stop, opt->verbosity, &eta, opt->threads);
			}
			opt->reset(opt->data);
			if(eta_adapt_result == OPT_SUCCESS){
				if(opt->verbosity > 0) printf("Stochastic gradient ascent using eta: %f\n", eta);
				if (opt->algorithm == OPT_SG_ADAM) {
					result = optimize_stochastic_gradient_adam(opt->maximize, opt->parameters, opt->f, opt->grad_f, eta, opt->data, &opt->stop, opt->verbosity, fmin, &checkpointer, opt->logger);
				}
				else{
					result = optimize_stochastic_gradient(opt->parameters, opt->f, opt->grad_f, eta, opt->data, &opt->stop, opt->verbosity, fmin, &checkpointer);
				}
				opt->reset(opt->data);
			}
			else{
				result = OPT_FAIL;
			}
			if (result == OPT_FAIL) {
				fprintf(stderr, "Stochastic gradient ascent failed\n");
			}
            break;
        }
		case OPT_TOPOLOGY:{
			result = topology_optimize(opt->data, fmin);
			break;
		}
		case OPT_EM:{
			result = em_optimize(opt, fmin);
			break;
		}
		default:
			result = -100;
			break;
	}
	if(opt->stop.oldx != NULL){
		free(opt->stop.oldx);
		opt->stop.oldx = NULL;  // reused across calls (bootstrap replicates, meta sweeps)
	}
	return result;
}

opt_result opt_optimize_univariate( Optimizer *opt, Parameter *p, double *fmin ){
	if ( opt->stop.time_max != 0 ) {
		time( &opt->stop.time_start );
	}
	opt_result result = OPT_SUCCESS;
	
	opt->dimension = 1;
	
	// should probably moved somewhere else
	opt->stop.oldx = dvector(opt->dimension);
    opt->stop.iter = 0;
    opt->stop.f_eval_current = 0;
    opt->stop.count = 0;
    opt->stop.stall = 0;
    
    if( opt->algorithm != OPT_BRENT ){
        error("optimize_univariate only works with Brent algorithm\n");
    }
    Parameters *ps = new_Parameters(1);
    Parameters_add(ps, p);
    
	result = brent_optimize( ps, opt->f, opt->data, &opt->stop, fmin );
	
    free_Parameters(ps);
	if(opt->stop.oldx != NULL){
		free(opt->stop.oldx);
		opt->stop.oldx = NULL;  // reused across calls (bootstrap replicates, meta sweeps)
	}
	return result;
}

opt_result opt_maximize( Optimizer *opt, double *fmin ){
	opt_result result = opt_optimize(opt, fmin);
    *fmin = - *fmin;
	return result;
}

opt_result opt_maximize_univariate( Optimizer *opt, Parameter *p, double *fmin ){
	opt_result result = opt_optimize_univariate(opt, p, fmin);
    *fmin = - *fmin;
	return result;
}

void opt_set_update_data_function( Optimizer *opt, opt_update_data uf ){
	if( opt != NULL ){
		opt->update = uf;
	}
}



bool xStop( const Parameters *x, double *xold, const double tolx){
	bool stop = true;
	
	for (int i = 0; i < Parameters_count(x); i++){
		//fprintf(stderr, "x=%f xold=%f tolx=%f (%d)\n", Parameters_value(x,i), xold[i], tolx,(fabs(Parameters_value(x,i) - xold[i] ) > tolx));
		if ( fabs(Parameters_value(x,i) - xold[i] ) > tolx ){
			stop = false;
		}
		xold[i] = Parameters_value(x,i);
	}
	
	return stop;
}



// algorithms that call opt->grad_f and therefore differentiate the target model
// with respect to opt->parameters. META is not one of them: each of its children
// is built by its own new_Optimizer_from_json call and checked there.
static bool opt_algorithm_uses_gradient(opt_algorithm algorithm){
	return algorithm == OPT_BFGS || algorithm == OPT_CG_PR || algorithm == OPT_CG_FR ||
	       algorithm == OPT_SG || algorithm == OPT_SG_ADAM;
}

Optimizer* new_Optimizer_from_json(json_node* node, Hashtable* hash){

	const char* algorithm_string = get_json_node_value_string(node, "algorithm");
	if (algorithm_string == NULL) {
		fprintf(stderr, "The `algorithm' key is not specified for object %s\n", get_json_node_value_string(node, "id"));
		exit(13);
	}
	
	if(strcasecmp(algorithm_string, "topology") == 0){
		Optimizer* opt = new_Optimizer(OPT_TOPOLOGY);
		opt->data = new_TopologyOptimizer_from_json(node, hash);
		return opt;
	}
	
	static const json_field schema[] = {
	    {"algorithm", JSON_REQUIRED, JSON_STRING},
	    {"alpha", JSON_OPTIONAL, JSON_NUMBER},
	    {"checkpoint", JSON_OPTIONAL, JSON_ANY},
	    {"checkpoint_frequency", JSON_OPTIONAL, JSON_NUMBER},
	    {"eta", JSON_OPTIONAL, JSON_NUMBER},
	    {"evaluations", JSON_OPTIONAL, JSON_NUMBER},
	    {"frequency_check", JSON_OPTIONAL, JSON_NUMBER},
	    {"iterations", JSON_OPTIONAL, JSON_NUMBER},
	    {"list", JSON_OPTIONAL, JSON_ANY},
	    {"logger", JSON_OPTIONAL, JSON_OBJECT},
	    {"max", JSON_FORBIDDEN, JSON_ANY},
	    {"maximize", JSON_OPTIONAL, JSON_ANY},
	    {"min", JSON_OPTIONAL, JSON_ANY},
	    {"model", JSON_FORBIDDEN, JSON_STRING|JSON_OBJECT},
	    {"parameters", JSON_OPTIONAL, JSON_ANY},
	    {"patience", JSON_OPTIONAL, JSON_NUMBER},
	    {"precision", JSON_OPTIONAL, JSON_NUMBER},
	    {"rounds", JSON_OPTIONAL, JSON_ANY},
		{"target", JSON_REQUIRED, JSON_STRING|JSON_OBJECT},
	    {"threads", JSON_OPTIONAL, JSON_NUMBER},
	    {"time", JSON_OPTIONAL, JSON_NUMBER},
	    {"tol", JSON_OPTIONAL, JSON_NUMBER},
	    {"treelikelihood", JSON_OPTIONAL, JSON_STRING},
	    {"update", JSON_OPTIONAL, JSON_STRING},
	    {"verbosity", JSON_OPTIONAL, JSON_NUMBER},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
	json_validate_xor(node, "treelikelihood", "parameters", "list", NULL);

	const char* id = get_json_node_value_string(node, "id");
	size_t iterations = get_json_node_value_size_t(node, "iterations", 1000);
	size_t min = get_json_node_value_size_t(node, "min", 1);
	bool maximize = get_json_node_value_bool(node, "maximize", true);
	double precision = get_json_node_value_double(node, "precision", 0.001);
	Parameters* parameters = new_Parameters(1);
	Optimizer* opt = NULL;
	json_node* parametersNode = get_json_node(node, "parameters");

	if(get_json_node(node, "treelikelihood") != NULL && parametersNode != NULL){
		fprintf(stderr, "Cannot specify both `treelikelihood' and `parameters' for object %s\n", id);
		exit(13);
	}
	if(parametersNode != NULL){
		get_parameters_references(node, hash, parameters);
	}
	
	if (strcasecmp(algorithm_string, "meta") == 0) {
		opt = new_Optimizer(OPT_META);
		OptimizerSchedule* schedule = opt_get_schedule(opt);
		json_node* list_node = get_json_node(node, "list");
		for (int i = 0; i < list_node->child_count; i++) {
			json_node* child = list_node->children[i];
			Optimizer* opt_child = new_Optimizer_from_json(child, hash);
			opt_add_optimizer(opt, opt_child);
			int child_rounds = get_json_node_value_int(child, "rounds", 1);
			opt->schedule->rounds[opt->schedule->count-1] = child_rounds;
		}
		opt_set_tolfx(opt, precision);
	}
	// BFGS
	else if (strcasecmp(algorithm_string, "bfgs") == 0) {
		opt = new_Optimizer(OPT_BFGS);
		opt_set_tolfx(opt, precision);
		opt->etas = dvector(1);
		// Quasi-Newton: try the full Newton step (alpha = 1) first.
		opt->etas[0] = get_json_node_value_double(node, "alpha", 1.0);
		// clone_Optimizer copies etas by eta_count; leaving it at 0 handed the
		// clone a NULL etas that dfpmin_optimize then dereferences.
		opt->eta_count = 1;
	}
	// Conjugate gradient
	else if (strcasecmp(algorithm_string, "cg") == 0) {
		opt = new_Optimizer(OPT_CG_PR);
		opt_set_tolfx(opt, precision);
	}
	else if (strcasecmp(algorithm_string, "brent") == 0 || strcasecmp(algorithm_string, "serial") == 0) {
		json_node* treelike_node = get_json_node(node, "treelikelihood");
		if (treelike_node != NULL) {
			const char* ref = (char*)treelike_node->value;
			opt = new_Optimizer(OPT_SERIAL_BRENT);
			opt->treelikelihood = safe_get_reference_model(ref, hash, id);
			//opt->treelikelihood->ref_count++;
		}
		else{
			if (Parameters_count(parameters) == 1 && Parameter_size(Parameters_at(parameters, 0)) == 1){
				opt = new_Optimizer(OPT_BRENT);
			}
			else{
				opt = new_Optimizer(OPT_SERIAL_BRENT);
			}
		}
		opt_set_tolx(opt, precision);
	}
	// EM for a free-rate site model. "treelikelihood" is not optional here as it
	// is for serial Brent: it is where the mixture being split lives, and there
	// is no parameter list to fall back on -- EM decides for itself which
	// parameters it moves.
	else if (strcasecmp(algorithm_string, "em") == 0) {
		json_node* treelike_node = get_json_node(node, "treelikelihood");
		if (treelike_node == NULL) {
			json_die(node, "\"algorithm\": \"em\" needs the \"treelikelihood\" whose "
			               "site model holds the mixture to split");
		}
		opt = new_Optimizer(OPT_EM);
		opt->treelikelihood = safe_get_reference_model((char*)treelike_node->value,
		                                               hash, id);
		opt_set_tolx(opt, precision);
		// "iterations" defaults to 1000, which is meaningless here: one EM step
		// already contains a full Brent run per category. Left at 0, EM applies
		// its own default of one step per category (IQ-TREE's).
		if (get_json_node(node, "iterations") == NULL) iterations = 0;
	}
    // stochastic gradient
    else if(strcasecmp(algorithm_string, "sg") == 0){
		char* update = get_json_node_value_string(node, "update");
		if (update != NULL) {
			if (strcasecmp(update, "adam") == 0) {
				opt = new_Optimizer(OPT_SG_ADAM);
				printf("Using adam\n");
			}
		}
		
		if(opt == NULL){
        	opt = new_Optimizer(OPT_SG);
		}
		opt->stop.frequency_check = get_json_node_value_size_t(node, "frequency_check", 100);
        opt->stop.tolfx = get_json_node_value_double(node, "tol", 0.001);
        json_node* etas = get_json_node(node, "eta");
        if (etas != NULL && etas->node_type == MJSON_ARRAY) {
			opt->etas = dvector(etas->child_count);
			for (int i = 0; i < etas->child_count; i++) {
				json_node* child = etas->children[i];
				opt->etas[i] = atof((char*)child->value);
			}
			opt->eta_count = etas->child_count;
		}
		else{
			opt->etas = dvector(1);
			opt->etas[0] = get_json_node_value_double(node, "eta", 1.0);
			opt->eta_count = 1;
		}
    }
	
	// Every branch above may leave opt NULL: an algorithm string that matches
	// none of them used to fall straight through to the dereference below.
	if (opt == NULL) {
		fprintf(stderr, "%s - unknown optimizer algorithm `%s'\n", id, algorithm_string);
		exit(13);
	}

	opt->maximize = maximize;
	opt_set_max_iteration(opt, iterations);
	opt_set_min_iteration(opt, min);
	opt->stop.patience = get_json_node_value_size_t(node, "patience", 1);

	// Resource budgets. Both default to 0, meaning no limit. The meta optimizer
	// enforces them between sweeps of its schedule, so they are honoured to
	// within one sweep; "evaluations" only counts what the schedule's entries
	// report (see opt_check_limits).
	double time_max = get_json_node_value_double(node, "time", 0);
	if (time_max > 0) opt_set_time_max(opt, time_max);
	size_t max_evaluations = get_json_node_value_size_t(node, "evaluations", 0);
	if (max_evaluations > 0) opt_set_max_evaluation(opt, max_evaluations);

	json_node* target_node = get_json_node(node, "target");
	Model* model = NULL;
	if(target_node->node_type == MJSON_STRING){
		const char* ref = (char*)target_node->value;
		model = safe_get_reference_model(ref, hash, id);
	}
	else{
		model = model_factory_from_json(target_node, hash);
	}
	opt_set_data(opt, model);
	
	if(maximize){
		opt_set_objective_function(opt, model_negative_logP);
		opt->grad_f = _negative_gradient;
	}
	else{
		opt_set_objective_function(opt, model_logP);
		opt->grad_f = _gradient;
	}
	opt->reset = _reset;
    
	if(parametersNode != NULL){
		if(opt_algorithm_uses_gradient(opt->algorithm)){
			Parameters_check_leaves(parameters, id);
		}
		opt_set_parameters(opt, parameters);
	}

	opt->threads = get_json_node_value_size_t(node, "threads", 1);
	opt->verbosity = get_json_node_value_int(node, "verbosity", 1);

	const char* checkpoint_file = get_json_node_value_string(node, "checkpoint");
	
	if (checkpoint_file != NULL){
		opt->checkpoint_file = String_clone(checkpoint_file);		
	}
	else{
		opt->checkpoint_file = String_clone("checkpoint.csv");
	}
	opt->checkpoint_frequency = get_json_node_value_int(node, "checkpoint_frequency", 1000);
	
	json_node* loggerNode = get_json_node(node, "logger");
	if(loggerNode != NULL){
		opt->logger = new_Trace_from_json(loggerNode, hash);
	}
	free_Parameters(parameters);
	return opt;
}

