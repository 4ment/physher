// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "cat.h"

#include <strings.h>

#include "parameters.h"
#include "treelikelihood.h"
#include "utils.h"
#include "matrix.h"


#define CAT_PROBE_DEFAULT 20
#define CAT_LLOYD_ITERATIONS 100
// Half the log-rate span the refined band of the probe grid extends past the
// previous round's centres on each side, as a fraction of that span.
#define CAT_WARM_MARGIN 0.5
// Floor on that margin, in log rate: log 2, so a round that left every pattern in
// one category still gets a band spanning a factor of four rather than a point.
#define CAT_WARM_MARGIN_MIN 0.6931471805599453
// Probe rates kept outside the band on each side, to hold the tails of the grid
// open at low resolution. Three is enough to keep a tail spanning a decade or two
// from being represented by its endpoint alone.
#define CAT_WARM_TAIL 3
// Fewest probe rates a refined grid is worth building on: below this the tails eat
// the budget and the band is no finer than the fixed grid it replaces.
#define CAT_WARM_MIN_PROBES 8
// Smallest rate the grid will probe at, RAxML's floor on its own per-pattern rate.
#define CAT_RATE_FLOOR 1.e-4
// Relative slack on the guard below. A reassignment that reproduces the model it
// started from -- every pattern in one category is the single-rate model, whichever
// category that is -- comes back a few ULPs apart because the partials are summed in
// a different order, and reverting on that would refuse a step that changed nothing.
// The same idiom, and the same order of magnitude, as the meta schedule's tolfx test.
#define CAT_REVERT_TOLERANCE 1.e-10
// Passes the nonparametric prior is allowed, and the Kiefer-Wolfowitz gap it is
// happy to stop at. The pass never touches the tree, so both are cheap in absolute
// terms -- P x M arithmetic against a tree traversal's P x nodes x states^2 -- but
// EM crawls once the estimate is nearly sparse, and on fluA the gap is still around
// 5.e-5 after ten thousand passes. The cap is therefore what normally ends the
// iteration, and it is set where the answer has stopped moving rather than where
// the certificate is satisfied: on that alignment the tree length is stable to five
// figures by pass 500 and the CAT likelihood to a fifth of a log unit by pass 1000,
// while a further nine thousand passes buy 0.16 of one. A run that wants the
// certificate raises "npmle_iterations".
#define CAT_NPMLE_ITERATIONS 1000
#define CAT_NPMLE_TOLERANCE 1.e-4
// Floor on the mass the estimated prior may put on a grid rate. The nonparametric
// MLE genuinely wants zeros -- that is what makes it discrete -- but a zero read as
// a prior is an unbreakable veto, and one arrived at by an iteration that was cut
// off at CAT_NPMLE_ITERATIONS is a veto nobody chose. The floor bounds it instead:
// a pattern has to beat the supported rates by log(1/CAT_NPMLE_FLOOR) = 230 log
// units before a floored rate can win, which no single alignment column does. It
// also keeps the mixture density strictly positive, so the EM ratio below cannot
// divide by zero however sparse the estimate becomes.
#define CAT_NPMLE_FLOOR 1.e-100
// What counts as an atom when reporting the support: mass above this is a rate the
// estimate actually uses, rather than one it is on its way to discarding.
#define CAT_NPMLE_ATOM 1.e-6

CatOptions cat_options_default(cat_assignment_t assignment){
	CatOptions options;
	options.assignment = assignment;
	options.prior = CAT_PRIOR_FIXED;
	options.npmle_iterations = CAT_NPMLE_ITERATIONS;
	options.npmle_tolerance = CAT_NPMLE_TOLERANCE;
	// FastTree hard-codes a Gamma(3, 1/3) prior and the arg-max rule was tuned
	// against it, so that pairing is kept. The posterior mean already shrinks by
	// the right amount on its own, and stacking the prior on top of it shrinks
	// twice: it costs 10 % of the tree length on large trees, so the default
	// there is a flat prior over the probe grid.
	options.prior_shape = (assignment == CAT_ASSIGNMENT_ARGMAX ? 3.0 : 0.0);
	options.probe_count = (assignment == CAT_ASSIGNMENT_ARGMAX ? 0 : CAT_PROBE_DEFAULT);
	options.verbosity = 0;
	return options;
}

// log of a mean-one Gamma(shape, 1/shape) density, up to a constant. The constant
// is the same for every category and only differences are ever used.
static void _cat_log_prior(double* log_prior, const double* rates, int count,
                           double shape){
	for (int i = 0; i < count; i++) {
		log_prior[i] = (shape - 1.0)*log(rates[i]) - shape*rates[i];
	}
}

// How the estimated prior turned out, for CatResult and the verbose print.
typedef struct {
	int atoms;
	int passes;
	double gap;
	// L(g) at the estimate: the log-likelihood of the mixture over the grid, which
	// is what the EM below maximizes. Comparable with a mixture's
	// score in a way the CAT score is not, and free -- the profile is already in
	// hand. Left at zero when the prior is not estimated.
	double logP;
} _cat_npmle_report;

// The per-pattern posterior over the grid the profile was built on, reduced to the
// numbers that say whether the grid is finer than the columns can resolve. All of
// it comes out of the array the assignment already built, so it costs no traversal.
// The fields are CatResult's; see cat.h for what each one means.
typedef struct {
	int grid;
	double confidence;
	double entropy;
	double information;
} _cat_diagnostics;

static void _cat_diagnose(const double* likelihoods, const double* log_prior,
                          int count, const double* pattern_weights,
                          int pattern_count, _cat_diagnostics* diagnostics){
	double sites = 0;
	for (int i = 0; i < pattern_count; i++) sites += pattern_weights[i];
	double* posterior = dvector(count);
	// Site-weighted aggregate posterior, whose entropy is the H(Z) half of the
	// mutual information below.
	double* aggregate = dvector(count);
	double confidence = 0;
	double entropy = 0;
	for (int i = 0; i < pattern_count; i++) {
		double max = -DBL_MAX;
		for (int j = 0; j < count; j++) {
			double value = likelihoods[j*pattern_count + i] + log_prior[j];
			if (value > max) max = value;
		}
		double total = 0;
		for (int j = 0; j < count; j++) {
			posterior[j] = exp(likelihoods[j*pattern_count + i] + log_prior[j] - max);
			total += posterior[j];
		}
		double best = 0;
		double row = 0;
		for (int j = 0; j < count; j++) {
			double p = posterior[j]/total;
			if (p > best) best = p;
			if (p > 0) row -= p*log(p);
			aggregate[j] += pattern_weights[i]*p;
		}
		confidence += pattern_weights[i]*best;
		entropy += pattern_weights[i]*row;
	}
	diagnostics->grid = count;
	diagnostics->confidence = (sites > 0 ? confidence/sites : 0);
	diagnostics->entropy = (sites > 0 ? entropy/sites : 0);
	double aggregate_entropy = 0;
	for (int j = 0; j < count; j++) {
		double p = (sites > 0 ? aggregate[j]/sites : 0);
		if (p > 0) aggregate_entropy -= p*log(p);
	}
	// H(Z) - H(Z | D) with both taken on the same grid: non-negative in exact
	// arithmetic, and clamped so a cancellation cannot report a negative
	// information at a grid nothing distinguishes.
	diagnostics->information = fmax(aggregate_entropy - diagnostics->entropy, 0.0);
	free(posterior);
	free(aggregate);
}

static void _cat_report_diagnostics(const _cat_diagnostics* diagnostics){
	printf("CAT assignment over %d grid rates: mean posterior %g, mean entropy %g "
	       "of %g nats, information %g\n", diagnostics->grid,
	       diagnostics->confidence, diagnostics->entropy, log(diagnostics->grid),
	       diagnostics->information);
}

// Nonparametric maximum likelihood estimate of the prior over the rate grid.
//
// Both assignment rules weight the per-pattern profile by a prior over the grid,
// and until now that prior was picked rather than estimated: flat, or FastTree's
// mean-one Gamma(3, 1/3). Neither is the rate distribution of the alignment in
// front of it, and the posterior mean is only an empirical-Bayes estimator if the
// prior it shrinks toward came from the data. Estimating it is the g-modelling half
// of empirical Bayes, and it can be done without assuming a family at all: over
// distributions g on the M grid rates the marginal log-likelihood
//
//     L(g) = sum_i n_i log( sum_j g_j L_ij )
//
// is concave, so it has a well-defined maximiser -- the classical nonparametric MLE
// of a mixing distribution (Kiefer & Wolfowitz 1956; Laird 1978; Koenker & Mizera
// 2014). L_ij is the profile the assignment already built, so the whole estimate
// costs no tree traversals, only M x P arithmetic per pass.
//
// The pass is the EM (vertex-direction) update
//
//     D_j = sum_i n_i L_ij / sum_k g_k L_ik,      g_j <- g_j D_j / N,
//
// which stays on the simplex -- sum_j g_j D_j / N = (1/N) sum_i n_i = 1 -- and
// cannot decrease L. D_j is also the certificate: the Kiefer--Wolfowitz condition
// for optimality is D_j <= N at every grid rate, with equality wherever g puts mass,
// so sup_j D_j/N - 1 measures how far the current g is from the maximum rather than
// how far the last step moved. That gap is what the loop stops on, and what it
// reports when it runs out of passes instead.
//
// Two properties of the answer are worth expecting rather than being surprised by.
//
// It is *discrete* -- the maximiser is supported on at most P of the M rates, and
// the mass on the rest goes to zero -- but it is approached from a dense start and
// the doomed masses only decay geometrically, so `report->atoms` after a bounded
// number of passes measures how far the iteration got as much as it measures the
// true support. On fluA it is 18 of 20 after a thousand passes. What is diagnostic
// is a *collapse*: an estimate that puts everything on one rate, and in particular
// on an end of the grid, is saying that no mixture beats rescaling the whole tree,
// which is what a badly scaled set of branch lengths looks like from here.
//
// And it is a prior on the rates, not a set of category weights: CAT assigns every
// pattern to one category with proportion 1, so what the estimate changes is where
// each pattern's posterior mean is pulled to, not how the likelihood mixes. Taking
// the atoms themselves as the categories -- the other half of what a nonparametric
// estimate offers -- would replace the quantizer, and is not what this does.
static void _cat_npmle(double* log_prior, const double* likelihoods,
                       const double* pattern_weights, int pattern_count,
                       const double* widths, int count, int max_passes,
                       double tolerance, _cat_npmle_report* report){
	// exp(profile - per-pattern maximum), so a pattern sitting at -800 does not
	// underflow its whole row to zero. The shift is constant along a row, which
	// scales that row's mixture density by a constant and so moves L(g) by an
	// additive constant: it cancels out of D_j exactly. The row's maximum entry is
	// 1, which is what keeps the mixture density strictly positive below.
	double* scaled = dvector((size_t)pattern_count*count);
	// The shift each row was scaled by, kept so L(g) can be put back on the scale
	// of a log-likelihood at the end.
	double* shift = dvector(pattern_count);
	for (int i = 0; i < pattern_count; i++) {
		double max = -DBL_MAX;
		for (int j = 0; j < count; j++) {
			double value = likelihoods[j*pattern_count + i];
			if (value > max) max = value;
		}
		shift[i] = max;
		for (int j = 0; j < count; j++) {
			scaled[(size_t)i*count + j] = exp(likelihoods[j*pattern_count + i] - max);
		}
	}

	// Start from the prior the grid itself implies: flat in rate density, so each
	// point carries its cell width, which is 1 everywhere on an evenly spaced grid.
	// The problem is concave, so the starting point decides how long the iteration
	// takes and not where it goes.
	double* g = dvector(count);
	double* gradient = dvector(count);
	double total = 0;
	for (int j = 0; j < count; j++) total += widths[j];
	for (int j = 0; j < count; j++) g[j] = widths[j]/total;

	double sites = 0;
	for (int i = 0; i < pattern_count; i++) sites += pattern_weights[i];

	// The gap is measured on the g that is about to be updated and the loop leaves
	// as soon as it is small enough, so `gap` on the way out always belongs to the
	// g that is returned and `pass` is the number of updates that produced it.
	//
	// Note that the gap is a certificate and not a progress meter: L(g) rises every
	// pass, but sup_j D_j can rise too, and on fluA the gap at pass 1000 is
	// sometimes larger than at pass 500. What that says is that a run is nearly
	// always stopped by `max_passes`, which is why the default sits where the
	// estimate has settled rather than where the gap has.
	double gap = INFINITY;
	int pass = 0;
	while (true) {
		for (int j = 0; j < count; j++) gradient[j] = 0;
		for (int i = 0; i < pattern_count; i++) {
			const double* row = scaled + (size_t)i*count;
			double mixture = 0;
			for (int j = 0; j < count; j++) mixture += g[j]*row[j];
			double weight = pattern_weights[i]/mixture;
			for (int j = 0; j < count; j++) gradient[j] += weight*row[j];
		}
		gap = 0;
		for (int j = 0; j < count; j++) {
			double excess = gradient[j]/sites - 1.0;
			if (excess > gap) gap = excess;
		}
		if (gap <= tolerance || pass >= max_passes) break;
		for (int j = 0; j < count; j++) {
			g[j] = fmax(g[j]*gradient[j]/sites, CAT_NPMLE_FLOOR);
		}
		pass++;
	}

	report->atoms = 0;
	report->passes = pass;
	report->gap = gap;
	for (int j = 0; j < count; j++) {
		if (g[j] > CAT_NPMLE_ATOM) report->atoms++;
		log_prior[j] = log(g[j]);
	}

	// L(g) itself, undoing the per-row shift. This is a genuine likelihood of the
	// data -- the label is summed out, not selected -- under the mixture that puts
	// mass g on the grid rates, at the tree the profile was taken on. It is
	// what the iteration maximizes, and reporting it costs one more pass over an
	// array that is about to be freed.
	report->logP = 0;
	for (int i = 0; i < pattern_count; i++) {
		const double* row = scaled + (size_t)i*count;
		double mixture = 0;
		for (int j = 0; j < count; j++) mixture += g[j]*row[j];
		report->logP += pattern_weights[i]*(shift[i] + log(mixture));
	}

	free(scaled);
	free(shift);
	free(g);
	free(gradient);
}

// The log prior *mass* each grid rate carries, which is what both assignment rules
// weight the profile by. `widths` is the cell width of each rate in log rate, so a
// fixed prior stays a prior on the rate rather than on the grid points however
// unevenly they are laid out; it is 1 everywhere on an evenly spaced grid, where
// this reduces to the log density.
static void _cat_prior(double* log_prior, const double* rates, const double* widths,
                       int count, const double* likelihoods,
                       const double* pattern_weights, int pattern_count,
                       const CatOptions* options, _cat_npmle_report* report){
	if (options->prior == CAT_PRIOR_NPMLE) {
		_cat_npmle(log_prior, likelihoods, pattern_weights, pattern_count, widths,
		           count, options->npmle_iterations, options->npmle_tolerance,
		           report);
		if (options->verbosity > 0) {
			printf("CAT estimated prior: %d of %d rates carry mass, "
			       "Kiefer-Wolfowitz gap %g after %d passes, mixture "
			       "log-likelihood %f\n", report->atoms, count, report->gap,
			       report->passes, report->logP);
		}
		return;
	}
	if (options->prior_shape > 0) {
		_cat_log_prior(log_prior, rates, count, options->prior_shape);
	}
	else {
		for (int j = 0; j < count; j++) log_prior[j] = 0;
	}
	for (int j = 0; j < count; j++) log_prior[j] += log(widths[j]);
}

// Fill likelihoods[count x P] with the per-pattern log-likelihood at each rate.
// Called with sm->cat_count forced to 1, so every pattern is scored against
// matrix block 0 and one traversal per rate is enough.
static void _cat_probe(SingleTreeLikelihood* tlk, const double* rates, int count,
                       double* likelihoods){
	SiteModel* sm = tlk->sm;
	for (int i = 0; i < count; i++) {
		sm->cat_rates[0] = rates[i];
		SingleTreeLikelihood_update_all_nodes(tlk);
		tlk->calculate(tlk);
		memcpy(likelihoods + (i*tlk->sp->count), tlk->pattern_lk,
		       sizeof(double)*tlk->sp->count);
	}
}

typedef struct {
	double value;
	double weight;
} _cat_point;

static int _cat_compare_points(const void* a, const void* b){
	double x = ((const _cat_point*)a)->value;
	double y = ((const _cat_point*)b)->value;
	return (x > y) - (x < y);
}

// Weighted Lloyd (k-means) in one dimension, seeded at the weighted quantiles.
// The rate estimates are already sorted along a single axis, so this converges in
// a handful of passes and its optimum is the least-squares K-point quantizer of
// the estimated rate distribution -- the empirical counterpart of the fixed
// 1/K..K grid it replaces.
static void _cat_quantize(const double* x, const double* weights, int n,
                          double* centres, int* labels, int count){
	_cat_point* sorted = malloc(sizeof(_cat_point)*n);
	for (int i = 0; i < n; i++) {
		sorted[i].value = x[i];
		sorted[i].weight = weights[i];
	}
	qsort(sorted, n, sizeof(_cat_point), _cat_compare_points);

	// Seed at the weighted quantiles of the estimates, so the categories start out
	// where the sites are. Several seeds can land on the same value -- identical
	// columns share a profile and therefore an estimate, and a real alignment puts
	// a large atom on one of them -- and the duplicates are nudged apart below
	// rather than spread over the support. Seeding one category per *distinct*
	// value instead is the obvious alternative and measures worse: it costs 6 % of
	// the tree length at 8-16 taxa, because it spends categories on a sparse tail
	// that carries almost no sites (notes/cat.md, "A note on choosing the quantizer").
	double total = 0;
	for (int i = 0; i < n; i++) total += sorted[i].weight;
	double accumulated = 0;
	int placed = 0;
	for (int i = 0; i < n && placed < count; i++) {
		accumulated += sorted[i].weight;
		if (accumulated >= total*(placed + 0.5)/count) centres[placed++] = sorted[i].value;
	}
	// Fewer usable quantiles than categories: the spares sit just above the last one
	// so the centres stay ordered, stay distinct enough for Lloyd to move them
	// apart, and no category is left uninitialised.
	while (placed < count) {
		centres[placed] = centres[placed - 1] + 1e-3;
		placed++;
	}
	free(sorted);

	for (int i = 0; i < n; i++) labels[i] = -1;
	for (int iteration = 0; iteration < CAT_LLOYD_ITERATIONS; iteration++) {
		bool moved = false;
		for (int i = 0; i < n; i++) {
			int best = 0;
			double best_distance = fabs(x[i] - centres[0]);
			for (int j = 1; j < count; j++) {
				double distance = fabs(x[i] - centres[j]);
				if (distance < best_distance) {
					best_distance = distance;
					best = j;
				}
			}
			if (labels[i] != best) {
				labels[i] = best;
				moved = true;
			}
		}
		if (!moved) break;
		for (int j = 0; j < count; j++) {
			double sum_weight = 0;
			double sum_value = 0;
			for (int i = 0; i < n; i++) {
				if (labels[i] == j) {
					sum_weight += weights[i];
					sum_value += weights[i]*x[i];
				}
			}
			// An empty cluster keeps its centre: dropping it would renumber the
			// categories and the rate vector is indexed by category.
			if (sum_weight > 0) centres[j] = sum_value/sum_weight;
		}
	}
}

// The grid the posterior mean is built on, and the weight each of its rates
// carries in that mean.
//
// The first call has nothing to go on and lays the rates out log-spaced over
// [1/probe, probe], every weight 1. Later calls know where the previous round put
// its centres and spend the same budget of rates unevenly: a band around those
// centres gets all but CAT_WARM_TAIL rates per side, so the resolution goes where
// the sites are, and the tails keep the grid open at their old endpoints. This is
// what makes successive rounds refine rather than reproduce round one's fixed
// point against updated branch lengths, which RAxML gets from a warm start at
// `patratStored[i]` on a step size that shrinks with the round.
//
// What the grid must *not* do is narrow onto the centres, which is the obvious
// reading of a warm start and is wrong here. RAxML's grid is a search device --
// the per-pattern rate it converges to does not depend on where the search
// started -- while this one is the support the posterior mean integrates over, so
// dropping its tails drops the model's ability to call a column slow or fast at
// all. Re-centring in that literal sense costs 75 log-likelihood units on
// examples/fluA/JC69-CAT20-ML.json. Hence: the range only ever grows, and the
// band's own edges move outward with the centres so a round whose rates have
// spread comes back with a wider band rather than a trapped one.
//
// Uneven spacing would change the estimate on its own -- the sum over the grid is
// a discrete prior, and crowding rates into the band would put prior mass there --
// so each rate carries its cell width in log rate as a weight. That leaves the
// posterior mean invariant to how the rates are distributed, and it is exactly 1
// everywhere on the evenly spaced grid, which is why a first call reproduces what
// it always did, bit for bit.
static void _cat_probe_grid(double* probe_rates, double* probe_weights,
                            int probe_count, const double* previous,
                            int previous_count){
	for (int i = 0; i < probe_count; i++) probe_weights[i] = 1.0;

	double lower = 1.0/probe_count;
	double upper = probe_count;
	double band_lower = lower;
	double band_upper = upper;
	if (previous != NULL && probe_count >= CAT_WARM_MIN_PROBES) {
		double lo = previous[0];
		double hi = previous[0];
		bool usable = true;
		for (int i = 0; i < previous_count; i++) {
			if (!isfinite(previous[i]) || previous[i] <= 0.0) usable = false;
			if (previous[i] < lo) lo = previous[i];
			if (previous[i] > hi) hi = previous[i];
		}
		if (usable) {
			double margin = fmax(CAT_WARM_MARGIN*(log(hi) - log(lo)),
			                     CAT_WARM_MARGIN_MIN);
			band_lower = fmax(exp(log(lo) - margin), CAT_RATE_FLOOR);
			band_upper = exp(log(hi) + margin);
			// The rates are normalised to a weighted mean of one, so a centre of
			// at least one is always among them and a band around them straddles
			// it; the tests are here so a grid can never come out of this
			// collapsed or inverted.
			if (!(band_lower < band_upper)) {
				band_lower = lower;
				band_upper = upper;
			}
			lower = fmin(lower, band_lower);
			upper = fmax(upper, band_upper);
		}
	}

	int tail_lower = (band_lower > lower ? CAT_WARM_TAIL : 0);
	int tail_upper = (band_upper < upper ? CAT_WARM_TAIL : 0);
	int band_count = probe_count - tail_lower - tail_upper;
	if (tail_lower == 0 && tail_upper == 0) {
		log_spaced_spaced_vector2(probe_rates, lower, upper, probe_count);
		return;
	}

	// The tails are laid out over the whole of their side and the endpoint the band
	// already covers is dropped, so no rate is duplicated and the outermost rate is
	// still the endpoint of the range.
	if (tail_lower > 0) {
		double* tail = dvector(tail_lower + 1);
		log_spaced_spaced_vector2(tail, lower, band_lower, tail_lower + 1);
		memcpy(probe_rates, tail, sizeof(double)*tail_lower);
		free(tail);
	}
	log_spaced_spaced_vector2(probe_rates + tail_lower, band_lower, band_upper,
	                          band_count);
	if (tail_upper > 0) {
		double* tail = dvector(tail_upper + 1);
		log_spaced_spaced_vector2(tail, band_upper, upper, tail_upper + 1);
		memcpy(probe_rates + tail_lower + band_count, tail + 1,
		       sizeof(double)*tail_upper);
		free(tail);
	}

	// Cell width in log rate, the outermost cells extended to the full width of the
	// one beside them rather than halved, so an evenly spaced grid comes out with
	// every weight equal. Scaled to a maximum of one: only ratios of weights are
	// ever used, and keeping them O(1) keeps them out of the way of the posterior's
	// own dynamic range.
	double largest = 0;
	for (int i = 0; i < probe_count; i++) {
		int left = (i == 0 ? 0 : i - 1);
		int right = (i == probe_count - 1 ? probe_count - 1 : i + 1);
		double width = (log(probe_rates[right]) - log(probe_rates[left]))
		               /(right - left);
		probe_weights[i] = width;
		if (width > largest) largest = width;
	}
	if (largest > 0) {
		for (int i = 0; i < probe_count; i++) probe_weights[i] /= largest;
	}
}

// Posterior-mean rate per pattern, then K categories at the cluster centres.
// Both halves are standard: the per-site posterior mean is the empirical-Bayes
// rate estimator (Mayrose et al. 2004), and clustering it is the least-squares
// quantizer (Lloyd 1982). What they buy together is that the categories are
// placed where the data are instead of on a fixed grid, and that each pattern is
// shrunk toward the mean by exactly as much as its column is uninformative.
// Returns the number of traversals _cat_probe performed.
static int _cat_assign_posterior_mean(SingleTreeLikelihood* tlk, double* rates,
                                      int cat_count, const CatOptions* options,
                                      const double* previous, int previous_count,
                                      _cat_npmle_report* report,
                                      _cat_diagnostics* diagnostics){
	int pattern_count = tlk->sp->count;
	int probe_count = options->probe_count;
	if (probe_count < cat_count) probe_count = cat_count;

	// A wider grid than the categories will end up spanning: the estimate has to
	// be free to sit outside the range they occupy, which is what caps the
	// dispersion the grid can represent.
	double* probe_rates = dvector(probe_count);
	double* probe_weights = dvector(probe_count);
	_cat_probe_grid(probe_rates, probe_weights, probe_count, previous,
	                previous_count);
	if (options->verbosity > 0) {
		printf("CAT probe grid [%f, %f] over %d rates, %s\n", probe_rates[0],
		       probe_rates[probe_count - 1], probe_count,
		       previous == NULL ? "evenly spaced"
		                        : "refined around the previous round's centres");
	}
	double* likelihoods = malloc(sizeof(double)*probe_count*pattern_count);
	_cat_probe(tlk, probe_rates, probe_count, likelihoods);

	double* weights = dvector(pattern_count);
	for (int i = 0; i < pattern_count; i++) weights[i] = tlk->sp->weights[i];

	// The prior comes after the probe because an estimated one is estimated *from*
	// the probe: the nonparametric MLE reads the whole P x M profile. A fixed prior
	// does not care when it is built.
	double* log_prior = dvector(probe_count);
	_cat_prior(log_prior, probe_rates, probe_weights, probe_count, likelihoods,
	           weights, pattern_count, options, report);

	// Taken on the probe grid rather than on the K categories: the posterior the
	// assignment reads is the one over the grid, and the categories it ends up
	// reporting are cluster centres of a point estimate, which carries no
	// posterior of its own.
	_cat_diagnose(likelihoods, log_prior, probe_count, weights, pattern_count,
	              diagnostics);
	if (options->verbosity > 0) _cat_report_diagnostics(diagnostics);

	// The quantizer works on log rates: the profile is close to symmetric there,
	// and the rate distributions this approximates (gamma, lognormal) put their
	// structure on the log scale.
	double* estimates = dvector(pattern_count);
	for (int i = 0; i < pattern_count; i++) {
		double max = -DBL_MAX;
		for (int j = 0; j < probe_count; j++) {
			double logP = likelihoods[j*pattern_count + i] + log_prior[j];
			if (logP > max) max = logP;
		}
		double numerator = 0;
		double denominator = 0;
		for (int j = 0; j < probe_count; j++) {
			// log_prior is the prior mass of the rate, cell width included, so a
			// grid whose rates are not evenly spaced still integrates the same
			// posterior.
			double posterior = exp(likelihoods[j*pattern_count + i] + log_prior[j] - max);
			numerator += posterior*probe_rates[j];
			denominator += posterior;
		}
		estimates[i] = log(numerator/denominator);
	}

	double* centres = dvector(cat_count);
	int* labels = ivector(pattern_count);
	_cat_quantize(estimates, weights, pattern_count, centres, labels, cat_count);
	for (int i = 0; i < cat_count; i++) rates[i] = exp(centres[i]);
	for (int i = 0; i < pattern_count; i++) tlk->sm->site_category[i] = labels[i];

	free(probe_rates);
	free(probe_weights);
	free(log_prior);
	free(likelihoods);
	free(estimates);
	free(weights);
	free(centres);
	free(labels);
	return probe_count;
}

// FastTree's rule: score every pattern at each of the K grid rates and keep the
// best, under the Gamma prior when there is one.
// Returns the number of traversals _cat_probe performed.
static int _cat_assign_argmax(SingleTreeLikelihood* tlk, const double* rates,
                              int cat_count, const CatOptions* options,
                              int* counts, _cat_npmle_report* report,
                              _cat_diagnostics* diagnostics){
	int pattern_count = tlk->sp->count;
	double* likelihoods = malloc(sizeof(double)*pattern_count*cat_count);
	_cat_probe(tlk, rates, cat_count, likelihoods);

	// This grid is evenly spaced in log rate, so every rate stands for the same
	// cell and the widths are all one; they are passed anyway so that both rules
	// build their prior through the same call.
	double* widths = dvector(cat_count);
	double* weights = dvector(pattern_count);
	for (int j = 0; j < cat_count; j++) widths[j] = 1.0;
	for (int i = 0; i < pattern_count; i++) weights[i] = tlk->sp->weights[i];

	double* log_prior = dvector(cat_count);
	_cat_prior(log_prior, rates, widths, cat_count, likelihoods, weights,
	           pattern_count, options, report);

	// The categories are the grid here, so max_j pi_ij is the posterior of the
	// category the loop below is about to select.
	_cat_diagnose(likelihoods, log_prior, cat_count, weights, pattern_count,
	              diagnostics);
	if (options->verbosity > 0) _cat_report_diagnostics(diagnostics);

	for (int i = 0; i < pattern_count; i++) {
		int best = 0;
		double bestLnl = -DBL_MAX;
		for (int j = 0; j < cat_count; j++) {
			double logP = likelihoods[j*pattern_count + i] + log_prior[j];
			if(logP > bestLnl){
				bestLnl = logP;
				best = j;
			}
		}
		tlk->sm->site_category[i] = best;
		if(options->verbosity > 0){
			printf("CAT pattern %d rate index %d rate %f [%f\n", i, best, rates[best],
			       bestLnl);
			counts[best]++;
		}
	}
	free(widths);
	free(weights);
	free(log_prior);
	free(likelihoods);
	return cat_count;
}

// Assign every pattern one rate category, by the rule in options->assignment.
//
// The probe scores every pattern against the block the pattern is currently
// assigned to, so it is only meaningful while every pattern sits in category 0:
// sm->cat_count is forced to 1 for the duration, and the transition-matrix refresh
// loop (treelikelihood.c, SingleTreeLikelihood_update_Q) therefore rewrites block 0
// alone. A pattern left in a category >= 1 by an earlier call would be scored against
// a block that never changes between probes, its profile would be flat, and the
// arg-max would hand it back to category 0 -- collapsing the whole assignment on the
// second call. Clearing the assignment first restores the precondition the first call
// gets for free, and is what makes repeated calls (an assign/optimize loop) possible.
CatResult cat_assign(SingleTreeLikelihood* tlk, const CatOptions* options){
	SiteModel* sm = tlk->sm;
	int cat_count = sm->cat_count;
	int pattern_count = tlk->sp->count;
	Parameter* cat_rates = Parameters_at(sm->rates, 0);
	CatResult result;
	result.reverted = false;
	_cat_npmle_report report = {0, 0, 0.0, 0.0};
	_cat_diagnostics diagnostics = {0, 0.0, 0.0, 0.0};

	// The number the new assignment has to beat, taken before site_category is
	// cleared: after that every pattern reads block 0 and the value would be the
	// single-rate model's, not the one the caller handed in. This is also the
	// traversal that brings the partials up to date for the probe, and it leaves
	// sm->need_update false so _cat_probe's writes into cat_rates[0] survive.
	result.logP = tlk->calculate(tlk);
	result.evaluations = 1;
	const double logP_before = result.logP;
	int* stored_categories = clone_ivector(sm->site_category, pattern_count);
	double* stored_rates = dvector(cat_count);
	for (int i = 0; i < cat_count; i++) {
		stored_rates[i] = Parameter_value_at(cat_rates, i);
	}
	const bool stored_scale = tlk->scale;

	// The centres the previous round left behind, to re-centre the probe grid on
	// (section "2. Warm start" of notes/cat.md), or NULL on the first call.
	//
	// They are read off sm->cat_rates rather than off the parameter because the
	// probe scores absolute rates: sm->cat_rates is the mean-one scale the branch
	// lengths were fitted against, the calculate above brought it up to date, and
	// _cat_probe is about to overwrite element 0 of it.
	//
	// Only the categories a pattern actually sits in are collected. An empty
	// cluster keeps a stale centre, and letting one pin the range would hold the
	// grid open for a rate no site uses. That doubles as the test for whether a
	// round has run at all: before the first call every pattern is in category 0,
	// and there is nothing to warm-start from.
	double* previous_centres = NULL;
	int previous_count = 0;
	int* occupied = ivector(cat_count);
	bool assigned = false;
	for (int i = 0; i < pattern_count; i++) occupied[sm->site_category[i]] = 1;
	for (int i = 1; i < cat_count; i++) assigned = assigned || occupied[i] != 0;
	if (assigned) {
		previous_centres = dvector(cat_count);
		for (int i = 0; i < cat_count; i++) {
			if (occupied[i] != 0) previous_centres[previous_count++] = sm->cat_rates[i];
		}
	}
	free(occupied);

	memset(sm->site_category, 0, sizeof(int)*pattern_count);
	sm->cat_count = 1;
	sm->need_update = false;
	int* counts = NULL;
	if(options->verbosity > 0) counts = ivector(cat_count);
	double* rates = dvector(cat_count);
	// The arg-max rule scores the categories themselves, so this grid *is* the
	// category set and the posterior mean overwrites it with its cluster centres.
	// It is rebuilt cold every call, and deliberately so for the arg-max: its
	// centres are its own grid endpoints by construction, so re-centring on them
	// the way the probe grid does would widen the span by a constant factor every
	// round instead of refining it. Warm-starting the arg-max needs RAxML's other
	// half -- a continuous per-pattern estimate to re-centre *on* -- which is
	// section 1, not section 2.
	log_spaced_spaced_vector2(rates, 1.0/cat_count, cat_count, cat_count);
//	log_spaced_spaced_vector2(rates+1, 1.0/(cat_count-1), cat_count-1, cat_count-1); // invariant
	if(options->verbosity > 0) print_dvector(rates, cat_count);

	if (options->assignment == CAT_ASSIGNMENT_POSTERIOR_MEAN) {
		result.evaluations += _cat_assign_posterior_mean(tlk, rates, cat_count, options,
		                                                previous_centres,
		                                                previous_count, &report,
		                                                &diagnostics);
		if(options->verbosity > 0){
			for (int i = 0; i < pattern_count; i++) counts[sm->site_category[i]]++;
			print_dvector(rates, cat_count);
		}
	}
	else {
		result.evaluations += _cat_assign_argmax(tlk, rates, cat_count, options, counts,
		                                         &report, &diagnostics);
	}

	// The category rates are one vector parameter. Every element but the last goes
	// in quietly and the last one loudly, so the listeners fire once, after the
	// whole vector is in place rather than on the first element of it.
	for (int i = 0; i < cat_count - 1; i++) {
		Parameter_set_value_at_quietly(cat_rates, rates[i], i);
	}
	Parameter_set_value_at(cat_rates, rates[cat_count - 1], cat_count - 1);
	if(options->verbosity > 0){
		for (int i = 0; i < cat_count; i++) {
			if(counts[i] == 0) fprintf(stdout, "CAT rate %f not used\n", rates[i]);
		}
		free(counts);
	}
	sm->cat_count = cat_count;
	SingleTreeLikelihood_use_rescaling(tlk, false);
	SingleTreeLikelihood_update_all_nodes(tlk);
	result.logP = tlk->calculate(tlk);
	result.evaluations++;
	free(rates);

	// Both halves of the step are heuristic -- the quantizer discards the spread
	// inside a category and the mean-one renormalization then moves every pattern
	// at once -- so neither is guaranteed to improve on the assignment that came
	// in. Put that one back when it does not, as RAxML's optimizeRateCategories
	// does (notes/cat.md, "4. The accept/revert guard"), so an alternation of this and the branch
	// lengths can only climb. Restoring the rescaling flag too is what makes the
	// revert exact: the accept path turns rescaling off, which is a different
	// numerical path through the same model.
	if (result.logP < logP_before - CAT_REVERT_TOLERANCE*fmax(1.0, fabs(logP_before))) {
		const double logP_proposed = result.logP;
		memcpy(sm->site_category, stored_categories, sizeof(int)*pattern_count);
		for (int i = 0; i < cat_count - 1; i++) {
			Parameter_set_value_at_quietly(cat_rates, stored_rates[i], i);
		}
		Parameter_set_value_at(cat_rates, stored_rates[cat_count - 1], cat_count - 1);
		SingleTreeLikelihood_use_rescaling(tlk, stored_scale);
		SingleTreeLikelihood_update_all_nodes(tlk);
		result.logP = tlk->calculate(tlk);
		result.evaluations++;
		result.reverted = true;
		if(options->verbosity > 0){
			fprintf(stdout, "CAT reassignment rejected: %f -> %f, restored %f\n",
			        logP_before, logP_proposed, result.logP);
		}
	}

	result.npmle_atoms = report.atoms;
	result.npmle_passes = report.passes;
	result.npmle_gap = report.gap;
	result.npmle_logP = report.logP;
	// Properties of the profile, so they describe the assignment that was proposed
	// whether or not the guard above kept it.
	result.grid = diagnostics.grid;
	result.confidence = diagnostics.confidence;
	result.entropy = diagnostics.entropy;
	result.information = diagnostics.information;

	free(stored_categories);
	free(stored_rates);
	free(previous_centres);
	return result;
}

// Sum the label out of the fit the model is currently holding.
//
// The CAT score is log P(D | zhat, theta) and every other site model reports
// log P(D | theta); the difference is the pointwise mutual information between a
// column and the category it was given, and no amount of tuning removes it
// (notes/cat.md, "What you may compare, and what you may not"). What does
// remove it is summing the label out, and the ingredients are already here: the
// category rates, and the share of sites each category holds.
//
// The mixture that comes out is not an approximation of CAT, it is CAT's own model
// read without the hard assignment -- a K-component free-rate model whose rates are
// the ones CAT fitted. Its mean rate is one for the same reason theirs is, because
// _cat_update normalised them against this very assignment, so the branch lengths
// it scores are the branch lengths in the tree and the number is comparable with a
// Gamma or free-rate score directly. What it is not comparable with is a Gamma
// score at equal parameter count: the mixture has 2K - 2 free rate parameters
// against Gamma's one, so an information criterion, not a raw difference, is the
// fair comparison.
//
// The profile is rebuilt here rather than borrowed from the assignment because the
// assignment probes the *raw* grid, before _cat_update divides through by the
// weighted mean; marginalizing that one would score a mixture at another tree
// scale. K traversals is the price of the honest version.
CatMixture cat_mixture(SingleTreeLikelihood* tlk, int verbosity){
	SiteModel* sm = tlk->sm;
	const int cat_count = sm->cat_count;
	const int pattern_count = tlk->sp->count;
	CatMixture mixture;
	mixture.evaluations = 0;

	// The CAT score, and the traversal that brings the partials up to date for the
	// probe. It leaves sm->need_update false, so _cat_probe's writes into
	// cat_rates[0] survive -- the same precondition cat_assign relies on.
	mixture.logP_cat = tlk->calculate(tlk);
	mixture.evaluations++;
	sm->update(sm);

	// The rates the likelihood just used, copied before the probe overwrites the
	// first of them. Read off sm->cat_rates rather than off the parameter because
	// the probe scores absolute rates and these are the mean-one ones; a site
	// model carrying a mu applies it to both alike, so it cancels.
	double* rates = clone_dvector(sm->cat_rates, cat_count);

	// Mixture weights: the share of sites each category holds. A hard assignment
	// has no soft responsibilities to average, and given the labels this is their
	// maximum-likelihood weight vector anyway.
	double* log_weight = dvector(cat_count);
	double sites = 0;
	for (int i = 0; i < pattern_count; i++) {
		log_weight[sm->site_category[i]] += tlk->sp->weights[i];
		sites += tlk->sp->weights[i];
	}
	mixture.used = 0;
	for (int k = 0; k < cat_count; k++) {
		if (log_weight[k] > 0) mixture.used++;
		// An empty category is a component of weight zero, and -inf is how it
		// drops out of the log-sum-exp below instead of distorting it.
		log_weight[k] = log(log_weight[k]/sites);
	}

	// _cat_probe's precondition: every pattern reading matrix block 0, and
	// sm->cat_count at 1 so the refresh loop rewrites that block alone.
	int* stored = clone_ivector(sm->site_category, pattern_count);
	memset(sm->site_category, 0, sizeof(int)*pattern_count);
	sm->cat_count = 1;
	sm->need_update = false;
	double* profile = dvector((size_t)cat_count*pattern_count);
	_cat_probe(tlk, rates, cat_count, profile);
	mixture.evaluations += cat_count;

	// Put back everything the probe borrowed. The category rates are rebuilt from
	// the parameter rather than from the copy above: same arithmetic on the same
	// inputs, so it lands on the same values, and it is the path any other caller
	// would take.
	sm->cat_count = cat_count;
	memcpy(sm->site_category, stored, sizeof(int)*pattern_count);
	sm->need_update = true;
	SingleTreeLikelihood_update_all_nodes(tlk);
	mixture.logP_cat = tlk->calculate(tlk);
	mixture.evaluations++;

	mixture.logP_mixture = 0;
	mixture.logP_bound = 0;
	for (int i = 0; i < pattern_count; i++) {
		double max = -DBL_MAX;
		for (int k = 0; k < cat_count; k++) {
			double value = profile[k*pattern_count + i] + log_weight[k];
			if (value > max) max = value;
		}
		double total = 0;
		for (int k = 0; k < cat_count; k++) {
			total += exp(profile[k*pattern_count + i] + log_weight[k] - max);
		}
		mixture.logP_mixture += tlk->sp->weights[i]*(max + log(total));
		// The same sum with the mixture replaced by a point mass on the selected
		// category: an ELBO, hence a lower bound, for any assignment at all.
		mixture.logP_bound += tlk->sp->weights[i]
		                       *(profile[stored[i]*pattern_count + i]
		                         + log_weight[stored[i]]);
	}
	mixture.gap = mixture.logP_cat - mixture.logP_mixture;

	// The mixture itself, so it can be rebuilt elsewhere. Empty categories are
	// skipped: they are components of weight zero and no other model has a way to
	// express one. Printed to enough digits that a "distribution": "discrete" site
	// model given these numbers reproduces logP_mixture rather than approaching it.
	if (verbosity > 0) {
		for (int k = 0; k < cat_count; k++) {
			if (log_weight[k] > -INFINITY) {
				printf("CAT component %d: rate %.10g weight %.10g\n", k, rates[k],
				       exp(log_weight[k]));
			}
		}
	}

	free(rates);
	free(log_weight);
	free(profile);
	free(stored);
	return mixture;
}

// The three numbers side by side, in the order they should be read: the one that
// may be reported, the one that may not, and the free bound that brackets the
// first from below.
//
// "at most" on the parameter count is load-bearing: 2K-2 counts every component as
// distinct and every value as fitted, and neither holds here. See the note on
// CatMixture::used.
static void _cat_report_mixture(const CatMixture* mixture){
	fprintf(stdout, "CAT mixture log-likelihood %f over %d categories "
	                "(at most %d free rate parameters)\n", mixture->logP_mixture,
	        mixture->used, 2*mixture->used - 2);
	// The score can land either side of the mixture: it is above it whenever the
	// assignment is the mixture's own arg-max, and a rule that selects on anything
	// else -- a prior, a posterior mean, a quantizer -- can put it below. Below is
	// the more damning reading of the two, since it says the hard assignment is
	// losing to the mixture built out of its own rates, so the direction is
	// spelled out rather than left to the sign of a number.
	fprintf(stdout, "CAT score %f, %s the mixture by %f; lower bound %f\n",
	        mixture->logP_cat, mixture->gap >= 0 ? "above" : "below",
	        fabs(mixture->gap), mixture->logP_bound);
}

void cat_check_sitemodel(json_node* node, const SingleTreeLikelihood* tlk){
	if (tlk->sm == NULL || tlk->sm->site_category == NULL) {
		json_die(node, "the tree likelihood does not carry an empirical CAT site "
		               "model: give its \"sitemodel\" a \"categories\" count and a "
		               "rate vector of that size, with no rate distribution");
	}
}

CatOptions cat_options_from_json(json_node* node){
	const char* assignment = get_json_node_value_string(node, "assignment");
	cat_assignment_t rule = CAT_ASSIGNMENT_ARGMAX;
	if (assignment != NULL) {
		if (strcasecmp(assignment, "argmax") == 0) {
			rule = CAT_ASSIGNMENT_ARGMAX;
		}
		else if (strcasecmp(assignment, "posterior_mean") == 0) {
			rule = CAT_ASSIGNMENT_POSTERIOR_MEAN;
		}
		else {
			json_die(node, "\"assignment\" is \"argmax\" or \"posterior_mean\", "
			               "not \"%s\"", assignment);
		}
	}

	CatOptions options = cat_options_default(rule);

	// "prior" is either a number -- the shape of the mean-one Gamma prior on the
	// category rate -- or the string "npmle", which estimates the prior instead of
	// naming one.
	//
	// The shape is FastTree's value and on by default for the arg-max: without it
	// the selected rate is unshrunk, the assigned rates are over-dispersed and the
	// branches stretch -- badly so on small trees. Set to 0 to select on the
	// likelihood alone. The posterior mean defaults to 0 for the opposite reason:
	// it shrinks on its own.
	json_node* prior_node = get_json_node(node, "prior");
	if (prior_node != NULL && prior_node->node_type == MJSON_STRING) {
		const char* prior = (const char*)prior_node->value;
		if (strcasecmp(prior, "npmle") == 0) {
			options.prior = CAT_PRIOR_NPMLE;
		}
		else {
			json_die(node, "\"prior\" is a Gamma shape or the string \"npmle\", "
			               "not \"%s\"", prior);
		}
	}
	else {
		options.prior_shape = get_json_node_value_double(node, "prior",
		                                                 options.prior_shape);
	}

	options.probe_count = get_json_node_value_int(node, "probe", options.probe_count);
	options.verbosity = get_json_node_value_int(node, "verbosity", 0);
	options.npmle_iterations = get_json_node_value_int(node, "npmle_iterations",
	                                                   options.npmle_iterations);
	options.npmle_tolerance = get_json_node_value_double(node, "npmle_tolerance",
	                                                     options.npmle_tolerance);
	if (rule == CAT_ASSIGNMENT_ARGMAX && get_json_node(node, "probe") != NULL) {
		json_die(node, "\"probe\" sets the size of the grid the posterior mean is "
		               "built on; \"assignment\": \"argmax\" scores the categories "
		               "themselves and has no separate grid");
	}
	if (options.probe_count < 0) {
		json_die(node, "\"probe\" must be positive, not %d", options.probe_count);
	}
	// The two NPMLE knobs only mean anything under the estimated prior, and a run
	// that set them and got FastTree's Gamma would be silently ignoring them.
	if (options.prior != CAT_PRIOR_NPMLE
	    && (get_json_node(node, "npmle_iterations") != NULL
	        || get_json_node(node, "npmle_tolerance") != NULL)) {
		json_die(node, "\"npmle_iterations\" and \"npmle_tolerance\" tune the "
		               "estimated prior; ask for it with \"prior\": \"npmle\"");
	}
	if (options.npmle_iterations < 1) {
		json_die(node, "\"npmle_iterations\" must be positive, not %d",
		         options.npmle_iterations);
	}
	if (!(options.npmle_tolerance >= 0)) {
		json_die(node, "\"npmle_tolerance\" must not be negative, not %g",
		         options.npmle_tolerance);
	}
	return options;
}

void cat_estimator_from_json(json_node* node, Hashtable* hash){
	static const json_field schema[] = {
	    {"assign", JSON_OPTIONAL, JSON_BOOL},
	    {"assignment", JSON_OPTIONAL, JSON_STRING},
	    {"id", JSON_OPTIONAL, JSON_STRING},
	    {"mixture", JSON_OPTIONAL, JSON_BOOL},
	    {"model", JSON_REQUIRED, JSON_STRING},
	    {"npmle_iterations", JSON_OPTIONAL, JSON_NUMBER},
	    {"npmle_tolerance", JSON_OPTIONAL, JSON_NUMBER},
	    {"prior", JSON_OPTIONAL, JSON_NUMBER | JSON_STRING},
	    {"probe", JSON_OPTIONAL, JSON_NUMBER},
	    {"type", JSON_OPTIONAL, JSON_STRING},
	    {"verbosity", JSON_OPTIONAL, JSON_NUMBER},
	};
	json_validate(node, schema, sizeof(schema)/sizeof(schema[0]));

	char* ref = get_json_node_value_string(node, "model");
	Model* mtlk = Hashtable_get(hash, ref+1);
	if (mtlk == NULL || mtlk->type != MODEL_TREELIKELIHOOD) {
		json_die(node, "\"model\" must reference a tree likelihood: '%s'", ref);
	}
	SingleTreeLikelihood *tlk = mtlk->obj;
	cat_check_sitemodel(node, tlk);

	// Two things this node can do, and they are usually wanted at different points
	// of a run. Reassigning belongs wherever the fit needs it; summing the label
	// out is a report and belongs after the *last* optimizer action, because the number
	// worth printing is the one at the tree that will be written out. Hence
	// "assign": false, which turns the node into a rescore that changes nothing --
	// the shape a run should use to report a comparable likelihood.
	bool assign = get_json_node_value_bool(node, "assign", true);
	bool mixture = get_json_node_value_bool(node, "mixture", !assign);
	if (!assign && !mixture) {
		json_die(node, "\"assign\": false with \"mixture\": false leaves this "
		               "node nothing to do");
	}
	if (!assign) {
		static const char* assignment_keys[] = {"assignment", "prior", "probe",
		                                        "npmle_iterations",
		                                        "npmle_tolerance"};
		for (size_t i = 0; i < sizeof(assignment_keys)/sizeof(assignment_keys[0]); i++) {
			if (get_json_node(node, assignment_keys[i]) != NULL) {
				json_die(node, "\"%s\" configures the assignment, which "
				               "\"assign\": false turns off", assignment_keys[i]);
			}
		}
	}

	if (assign) {
		CatOptions options = cat_options_from_json(node);
		CatResult result = cat_assign(tlk, &options);
		if (result.reverted) {
			fprintf(stdout, "CAT: the reassignment scored worse than the assignment "
			                "it started from (%f); it was refused and nothing "
			                "changed\n", result.logP);
		}
	}
	if (mixture) {
		CatMixture result = cat_mixture(tlk, get_json_node_value_int(node,
		                                                            "verbosity", 0));
		_cat_report_mixture(&result);
	}
}
