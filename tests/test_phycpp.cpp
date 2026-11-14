// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

// Tests of the C++ wrapper (libphycpp): every likelihood it builds is checked
// against physher's own JSON construction, and every gradient against central
// finite differences taken through the wrapper's SetParameters.

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "phycpp/physher.hpp"

extern "C" {
#include "phyc/ctmcscale.h"
#include "phyc/hashtable.h"
#include "phyc/mjson.h"
}

// minunit returns string literals as char*, which C++ rejects; same idea here.
#define mu_assert(test, message) \
    if (!(test)) {               \
        return message;          \
    }
#define mu_run_test(test)                  \
    do {                                   \
        const char *message = test();      \
        tests_run++;                       \
        if (message) return message;       \
    } while (0)

static int tests_run = 0;

namespace {

// The right child of the root ("D") is a tip and comes first in the taxon list,
// which is the case where physher's branch order and the wrapper's differ.
const char *kNewick = "(((A:0.1,B:0.2):0.1,C:0.3):0.05,D:0.2);";
const std::vector<std::string> kTaxa = {"D", "A", "B", "C"};
const std::vector<std::pair<std::string, std::string>> kAlignment = {
    {"A", "ACGTACGTAACCGGTTACGTAGCTAGCTAGGATCCA"},
    {"B", "ACGTACTTAACCGGTTACGTAGCTAGCAAGGATCCA"},
    {"C", "ACGAACGTAACCGCTTACGTTGCTAGCTAGGATGCA"},
    {"D", "TCGTACGTAACAGGTTACCTAGCTTGCTAGCATCCA"},
};
const std::vector<double> kDates = {0.0, 0.0, 0.0, 0.0};

std::string Sequences() {
    std::string s;
    for (const auto &[taxon, seq] : kAlignment) {
        if (!s.empty()) s += ",";
        s += "\"" + taxon + "\":\"" + seq + "\"";
    }
    return s;
}

// logP of a tree likelihood built from JSON, on the same data and tree.
double JsonLogP(const std::string &substitution, const std::string &sitemodel) {
    Hashtable *hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);
    std::string content =
        "{\"id\":\"tlk\",\"type\":\"treelikelihood\",\"tipstates\":false,"
        "\"sitepattern\":{\"id\":\"patterns\",\"type\":\"sitepattern\","
        "\"datatype\":\"nucleotide\",\"alignment\":{\"id\":\"seqs\","
        "\"type\":\"alignment\",\"sequences\":{" +
        Sequences() + "}}},\"sitemodel\":" + sitemodel +
        ",\"substitutionmodel\":" + substitution +
        ",\"tree\":{\"id\":\"tree\",\"type\":\"tree\",\"newick\":\"" + kNewick +
        "\",\"keep_branch_lengths\":true,\"branch_lengths\":{\"id\":\"bl\","
        "\"type\":\"parameter\",\"dimension\":5,\"x\":0.1}}}";
    json_node *json = create_json_tree(content.c_str());
    Model *model = new_TreeLikelihoodModel_from_json(json, hash);
    double logP = model->logP(model);
    model->free(model);
    free_Hashtable(hash);
    json_free_tree(json);
    return logP;
}

// Central finite differences of f wrt every element of x.
std::vector<double> FiniteDifferences(const std::function<double()> &f,
                                      ModelInterface *m, std::vector<double> x,
                                      double eps = 1.e-6) {
    std::vector<double> grad(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        double xi = x[i];
        x[i] = xi + eps;
        m->SetParameters(x.data());
        double up = f();
        x[i] = xi - eps;
        m->SetParameters(x.data());
        double down = f();
        x[i] = xi;
        grad[i] = (up - down) / (2 * eps);
    }
    m->SetParameters(x.data());
    return grad;
}

bool Close(double a, double b, double tol = 1.e-5) {
    return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b));
}

bool CloseAll(const double *a, const std::vector<double> &b, const char *what) {
    for (size_t i = 0; i < b.size(); i++) {
        if (!Close(a[i], b[i])) {
            fprintf(stderr, "%s[%zu]: analytic %.10f, finite differences %.10f\n", what,
                    i, a[i], b[i]);
            return false;
        }
    }
    return true;
}

const std::vector<double> kFrequencies = {0.3, 0.2, 0.2, 0.3};

}  // namespace

static const char *test_unrooted_jc69() {
    UnRootedTreeModelInterface tree(kNewick, kTaxa);
    JC69Interface jc69;
    ConstantSiteModelInterface sitemodel(std::nullopt);
    TreeLikelihoodInterface tlk(kAlignment, &tree, &jc69, &sitemodel, std::nullopt);

    double expected = JsonLogP(
        "{\"id\":\"sm\",\"type\":\"substitutionmodel\",\"kind\":\"jc69\","
        "\"datatype\":\"nucleotide\",\"frequencies\":{\"id\":\"freqs\","
        "\"type\":\"Simplex\",\"x\":[0.25,0.25,0.25,0.25]}}",
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\"}");
    mu_assert(Close(tlk.LogLikelihood(), expected, 1.e-10), "JC69 logP differs from JSON");

    // Get/Set round trip, in the wrapper's own order
    std::vector<double> bl(tree.parameterCount_);
    tree.GetParameters(bl.data());
    mu_assert(Close(bl[1], 0.1, 1.e-12) && Close(bl[3], 0.3, 1.e-12),
              "branch lengths of A and C are not at their taxon index");
    // D is the pinned right child of the root: its slot carries the length of the
    // unrooted branch between D and the rest, 0.05 + 0.2
    mu_assert(Close(bl[0], 0.25, 1.e-12), "unrooted root branch not in the slot of D");
    tree.SetParameters(bl.data());
    mu_assert(Close(tlk.LogLikelihood(), expected, 1.e-10), "Set/Get round trip");

    // gradient in the same order as SetParameters
    std::vector<double> grad(tree.parameterCount_);
    tlk.Gradient(grad.data(), {GradientFlags::TREE_MODEL_BRANCH_LENGTHS});
    auto fd = FiniteDifferences([&] { return tlk.LogLikelihood(); }, &tree, bl);
    mu_assert(CloseAll(grad.data(), fd, "branch"), "JC69 branch length gradient");
    return nullptr;
}

static const char *test_unrooted_gtr_gamma_invariant() {
    UnRootedTreeModelInterface tree(kNewick, kTaxa);
    GTRInterface gtr({1.5, 3.0, 0.8, 1.2, 4.0}, kFrequencies);
    GammaSiteModelInterface sitemodel(0.7, 4, 0.2, std::nullopt);
    TreeLikelihoodInterface tlk(kAlignment, &tree, &gtr, &sitemodel, std::nullopt);
    mu_assert(sitemodel.parameterCount_ == 2, "+G+I has a shape and a proportion");
    mu_assert(sitemodel.GetCategoryCount() == 5, "4 gamma categories + invariant");

    double expected = JsonLogP(
        "{\"id\":\"sm\",\"type\":\"substitutionmodel\",\"kind\":\"gtr\","
        "\"datatype\":\"nucleotide\",\"frequencies\":{\"id\":\"freqs\","
        "\"type\":\"Simplex\",\"x\":[0.3,0.2,0.2,0.3]},\"rates\":{"
        "\"ac\":{\"id\":\"ac\",\"type\":\"parameter\",\"x\":1.5,\"lower\":0},"
        "\"ag\":{\"id\":\"ag\",\"type\":\"parameter\",\"x\":3.0,\"lower\":0},"
        "\"at\":{\"id\":\"at\",\"type\":\"parameter\",\"x\":0.8,\"lower\":0},"
        "\"cg\":{\"id\":\"cg\",\"type\":\"parameter\",\"x\":1.2,\"lower\":0},"
        "\"ct\":{\"id\":\"ct\",\"type\":\"parameter\",\"x\":4.0,\"lower\":0}}}",
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\",\"distribution\":\"gamma\","
        "\"categories\":4,\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\","
        "\"x\":0.7,\"lower\":0},\"proportion_invariant\":0.2}");
    mu_assert(Close(tlk.LogLikelihood(), expected, 1.e-10),
              "GTR+G4+I logP differs from JSON");

    auto logP = [&] { return tlk.LogLikelihood(); };
    std::vector<double> grad(32);

    // site model: shape then proportion invariant
    std::vector<double> site(sitemodel.parameterCount_);
    sitemodel.GetParameters(site.data());
    mu_assert(Close(site[0], 0.7, 1.e-12) && Close(site[1], 0.2, 1.e-12),
              "site model GetParameters");
    tlk.Gradient(grad.data(),
                 {GradientFlags::SITE_MODEL_PARAMETER, GradientFlags::SITE_MODEL_PINV});
    mu_assert(CloseAll(grad.data(), FiniteDifferences(logP, &sitemodel, site), "site"),
              "shape/pinv gradient");

    // substitution model: 5 rates then 4 frequencies
    std::vector<double> subst = {1.5, 3.0, 0.8, 1.2, 4.0, 0.3, 0.2, 0.2, 0.3};
    tlk.Gradient(grad.data(), {GradientFlags::SUBSTITUTION_MODEL_RATES,
                               GradientFlags::SUBSTITUTION_MODEL_FREQUENCIES});
    mu_assert(CloseAll(grad.data(), FiniteDifferences(logP, &gtr, subst), "gtr"),
              "GTR rates/frequencies gradient");

    // everything at once: branch lengths, site model, substitution model
    std::vector<double> bl(tree.parameterCount_);
    tree.GetParameters(bl.data());
    tlk.Gradient(grad.data(), {GradientFlags::TREE_MODEL_BRANCH_LENGTHS,
                               GradientFlags::SITE_MODEL_PARAMETER,
                               GradientFlags::SITE_MODEL_PINV,
                               GradientFlags::SUBSTITUTION_MODEL_RATES});
    mu_assert(CloseAll(grad.data(), FiniteDifferences(logP, &tree, bl), "branch"),
              "branch lengths in a joint gradient");
    mu_assert(CloseAll(grad.data() + bl.size(),
                       FiniteDifferences(logP, &sitemodel, site), "site"),
              "site model offset in a joint gradient");
    std::vector<double> fdRates = FiniteDifferences(logP, &gtr, subst);
    fdRates.resize(5);
    mu_assert(CloseAll(grad.data() + bl.size() + site.size(), fdRates, "rates"),
              "substitution offset in a joint gradient");
    return nullptr;
}

static const char *test_invariant_hky() {
    UnRootedTreeModelInterface tree(kNewick, kTaxa);
    HKYInterface hky(2.5, kFrequencies);
    InvariantSiteModelInterface sitemodel(0.3, std::nullopt);
    TreeLikelihoodInterface tlk(kAlignment, &tree, &hky, &sitemodel, std::nullopt);
    mu_assert(sitemodel.parameterCount_ == 1, "+I has one parameter");

    double expected = JsonLogP(
        "{\"id\":\"sm\",\"type\":\"substitutionmodel\",\"kind\":\"hky\","
        "\"datatype\":\"nucleotide\",\"frequencies\":{\"id\":\"freqs\","
        "\"type\":\"Simplex\",\"x\":[0.3,0.2,0.2,0.3]},\"rates\":{\"kappa\":{"
        "\"id\":\"kappa\",\"type\":\"parameter\",\"x\":2.5,\"lower\":0}}}",
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\",\"distribution\":\"discrete\","
        "\"categories\":1,\"proportion_invariant\":0.3}");
    mu_assert(Close(tlk.LogLikelihood(), expected, 1.e-10), "HKY+I logP differs from JSON");

    double rates[2], proportions[2];
    sitemodel.GetRates(rates);
    sitemodel.GetProportions(proportions);
    mu_assert(rates[0] == 0 && Close(proportions[0], 0.3, 1.e-12), "+I categories");

    auto logP = [&] { return tlk.LogLikelihood(); };
    std::vector<double> grad(8);
    std::vector<double> pinv = {0.3};
    tlk.Gradient(grad.data(), {GradientFlags::SITE_MODEL_PINV});
    mu_assert(CloseAll(grad.data(), FiniteDifferences(logP, &sitemodel, pinv), "pinv"),
              "+I gradient");

    std::vector<double> subst = {2.5, 0.3, 0.2, 0.2, 0.3};
    tlk.Gradient(grad.data(), {GradientFlags::SUBSTITUTION_MODEL_RATES,
                               GradientFlags::SUBSTITUTION_MODEL_FREQUENCIES});
    std::vector<double> fd = FiniteDifferences(logP, &hky, subst);
    mu_assert(CloseAll(grad.data(), {fd[0]}, "kappa"), "HKY kappa gradient");
    // HKY's closed-form P(t) assumes the frequencies sum to 1, so only the
    // derivative along the simplex is defined: compare it along e_i - e_0.
    for (size_t i = 2; i < 5; i++) {
        const double eps = 1.e-6;
        std::vector<double> x = subst;
        x[i] += eps;
        x[1] -= eps;
        hky.SetParameters(x.data());
        double up = logP();
        x[i] -= 2 * eps;
        x[1] += 2 * eps;
        hky.SetParameters(x.data());
        double down = logP();
        double directional = grad[i] - grad[1];
        if (!Close(directional, (up - down) / (2 * eps))) {
            fprintf(stderr, "frequency direction %zu: analytic %.10f, fd %.10f\n", i,
                    directional, (up - down) / (2 * eps));
            return "HKY frequency gradient along the simplex";
        }
    }
    hky.SetParameters(subst.data());
    return nullptr;
}

// Setting the heights of a time tree must reach the likelihood: it used to keep
// the branch lengths cached from the first evaluation.
static const char *test_time_tree_heights() {
    TimeTreeModelInterface tree(kNewick, kTaxa, kDates);
    StrictClockModelInterface clock(1.0, &tree);
    JC69Interface jc69;
    ConstantSiteModelInterface sitemodel(std::nullopt);
    TreeLikelihoodInterface tlk(kAlignment, &tree, &jc69, &sitemodel, &clock);

    std::vector<double> heights(tree.parameterCount_);
    tree.GetParameters(heights.data());
    double before = tlk.LogLikelihood();
    std::vector<double> scaled = heights;
    for (double &h : scaled) h *= 1.5;
    tree.SetParameters(scaled.data());
    double after = tlk.LogLikelihood();
    mu_assert(!Close(before, after, 1.e-6), "new heights did not change the likelihood");
    double full = tlk.GetModel()->full_logP(tlk.GetModel());
    mu_assert(Close(after, full, 1.e-10), "likelihood stale after SetParameters");
    tree.SetParameters(heights.data());
    mu_assert(Close(tlk.LogLikelihood(), before, 1.e-10), "heights not restored");
    return nullptr;
}

static const char *test_reparameterized_time_tree() {
    ReparameterizedTimeTreeModelInterface tree(kNewick, kTaxa, kDates,
                                               TreeTransformFlags::RATIO);
    StrictClockModelInterface clock(0.8, &tree);
    JC69Interface jc69;
    ConstantSiteModelInterface sitemodel(std::nullopt);
    TreeLikelihoodInterface tlk(kAlignment, &tree, &jc69, &sitemodel, &clock);

    // n-2 ratios + root height
    std::vector<double> params(tree.parameterCount_);
    tree.GetParameters(params.data());
    std::vector<double> copy = params;
    tree.GetParameters(copy.data());
    mu_assert(copy == params, "GetParameters must not modify the tree");

    // same heights through a plain time tree give the same likelihood
    std::vector<double> heights(tree.parameterCount_);
    tree.GetNodeHeights(heights.data());
    TimeTreeModelInterface plain(kNewick, kTaxa, kDates);
    StrictClockModelInterface plainClock(0.8, &plain);
    JC69Interface plainJc69;
    ConstantSiteModelInterface plainSitemodel(std::nullopt);
    TreeLikelihoodInterface plainTlk(kAlignment, &plain, &plainJc69, &plainSitemodel,
                                     &plainClock);
    plain.SetParameters(heights.data());
    mu_assert(Close(tlk.LogLikelihood(), plainTlk.LogLikelihood(), 1.e-10),
              "reparameterized and plain time trees disagree");

    auto logP = [&] { return tlk.LogLikelihood(); };
    std::vector<double> grad(tree.parameterCount_ + 1);
    tlk.Gradient(grad.data(),
                 {GradientFlags::TREE_MODEL_RATIO, GradientFlags::BRANCH_MODEL_RATE});
    mu_assert(CloseAll(grad.data(), FiniteDifferences(logP, &tree, params), "ratios"),
              "ratio/root height gradient");
    std::vector<double> rate = {0.8};
    mu_assert(CloseAll(grad.data() + params.size(),
                       FiniteDifferences(logP, &clock, rate), "rate"),
              "clock rate gradient");
    return nullptr;
}

static const char *test_ctmc_scale() {
    ReparameterizedTimeTreeModelInterface tree(kNewick, kTaxa, kDates,
                                               TreeTransformFlags::RATIO);
    // the prior is a product over the rates
    CTMCScaleModelInterface both({0.5, 2.0}, &tree);
    CTMCScaleModelInterface first({0.5}, &tree);
    CTMCScaleModelInterface second({2.0}, &tree);
    mu_assert(Close(both.LogLikelihood(), first.LogLikelihood() + second.LogLikelihood(),
                    1.e-10),
              "CTMC scale must use every rate");

    double rate[1] = {0.7};
    first.SetParameters(rate);
    double got[1];
    first.GetParameters(got);
    mu_assert(got[0] == 0.7, "CTMC Get/SetParameters");
    return nullptr;
}

static const char *all_tests() {
    mu_run_test(test_unrooted_jc69);
    mu_run_test(test_unrooted_gtr_gamma_invariant);
    mu_run_test(test_invariant_hky);
    mu_run_test(test_time_tree_heights);
    mu_run_test(test_reparameterized_time_tree);
    mu_run_test(test_ctmc_scale);
    return nullptr;
}

int main() {
    const char *result = all_tests();
    if (result != nullptr) {
        printf("FAILED: %s\n", result);
    } else {
        printf("ALL TESTS PASSED\n");
    }
    printf("Tests run: %d\n", tests_run);
    return result != nullptr;
}
