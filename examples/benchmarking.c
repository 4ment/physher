#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "phyc/datatype.h"
#include "phyc/demographicmodels.h"
#include "phyc/gradient.h"
#include "phyc/parameters.h"
#include "phyc/sequenceio.h"
#include "phyc/sitepattern.h"
#include "phyc/tree.h"
#include "phyc/treeio.h"
#include "phyc/treelikelihood.h"
#include "phyc/treetransform.h"
#include "phyc/matrix.h"

double mseconds(struct timespec start, struct timespec end) {
    return (end.tv_sec - start.tv_sec) * 1000. +
           (end.tv_nsec - start.tv_nsec) / 1000000.;
}

char* json_jc69_strict =
    "{ \
	\"id\":\"treelikelihood\", \
	\"type\": \"treelikelihood\", \
	\"reparameterized\": true, \
	\"sse\":true, \
	\"tipstates\": false, \
	\"sitepattern\": \"&patterns\", \
	\"substitutionmodel\":{ \
		\"id\":\"sm\", \
		\"type\":\"substitutionmodel\", \
		\"model\":\"jc69\", \
		\"datatype\":\"nucleotide\", \
		\"frequencies\":{ \
			\"id\":\"freqs\", \
			\"type\":\"parameter\", \
			\"x\":[0.25,0.25,0.25,0.25] \
		} \
	}, \
	\"sitemodel\":{ \
		\"id\": \"sitemodel\", \
		\"type\": \"sitemodel\" \
	 \
	}, \
	\"tree\": \"&tree\", \
	\"branchmodel\":{ \
		\"id\": \"bm\", \
		\"type\": \"branchmodel\", \
		\"model\": \"strict\", \
		\"tree\": \"&tree\", \
		\"rate\": { \
			\"id\": \"rate\", \
            \"type\": \"parameter\", \
            \"x\": {\"id\": \"rate2\", \"type\": \"parameter\", \"x\":-6.907755278982137}, \
            \"lower\":0 \
		} \
	} \
}";

char* json_jc69_unrooted =
    "{ \
	\"id\":\"treelikelihood\", \
	\"type\": \"treelikelihood\", \
	\"sse\":true, \
	\"tipstates\": false, \
	\"sitepattern\": \"&patterns\", \
	\"substitutionmodel\":{ \
		\"id\":\"sm\", \
		\"type\":\"substitutionmodel\", \
		\"model\":\"jc69\", \
		\"datatype\":\"nucleotide\", \
		\"frequencies\":{ \
			\"id\":\"freqs\", \
			\"type\":\"parameter\", \
			\"x\":[0.25,0.25,0.25,0.25] \
		} \
	}, \
	\"sitemodel\":{ \
		\"id\": \"sitemodel\", \
		\"type\": \"sitemodel\" \
	}, \
	\"tree\": \"&tree\" \
}";

char* json_gtr_unrooted =
    "{ \
	\"id\":\"treelikelihood\", \
	\"type\": \"treelikelihood\", \
	\"sse\":true, \
	\"tipstates\": false, \
	\"sitepattern\": \"&patterns\", \
	\"substitutionmodel\":{ \
		\"id\":\"sm\", \
		\"type\":\"substitutionmodel\", \
		\"model\":\"gtr\", \
		\"datatype\":\"nucleotide\", \
		\"_rates\":{ \
			\"ac\":{\"id\":\"ac\", \"type\":\"parameter\", \"x\":{\"id\":\"ac2\", \"type\": \"parameter\", \"x\":0}, \"lower\":0, \"upper\":\"infinity\"}, \
			\"ag\":{\"id\":\"ag\", \"type\":\"parameter\", \"x\":{\"id\":\"ag2\", \"type\": \"parameter\", \"x\":0}, \"lower\":0, \"upper\":\"infinity\"}, \
			\"at\":{\"id\":\"at\", \"type\":\"parameter\", \"x\":{\"id\":\"at2\", \"type\": \"parameter\", \"x\":0}, \"lower\":0, \"upper\":\"infinity\"}, \
			\"cg\":{\"id\":\"cg\", \"type\":\"parameter\", \"x\":{\"id\":\"cg2\", \"type\": \"parameter\", \"x\":0}, \"lower\":0, \"upper\":\"infinity\"}, \
			\"ct\":{\"id\":\"ct\", \"type\":\"parameter\", \"x\":{\"id\":\"ct2\", \"type\": \"parameter\", \"x\":0}, \"lower\":0, \"upper\":\"infinity\"} \
		}, \
		\"rates\":{ \
			\"id\":\"rates\", \
			\"type\":\"Simplex\", \
			\"x\":{ \
                \"id\":\"rates2\", \
                \"type\": \"parameter\", \
                \"x\": [0.0, 0.0, 0.0, 0.0, 0.0] \
            } \
		}, \
		\"frequencies\":{ \
			\"id\":\"freqs\", \
			\"type\":\"Simplex\", \
			\"x\":{ \
                \"id\":\"freqs2\", \
                \"type\": \"parameter\", \
                \"x\": [0.0, 0.0, 0.0] \
            } \
		} \
	}, \
	\"sitemodel\":{ \
		\"id\": \"sitemodel\", \
		\"type\": \"sitemodel\" \
	 \
	}, \
	\"tree\": \"&tree\" \
}";

// Calculate the log of the determinant of the Jacobian for the node height
// transform. Derivatives are wrt ratio/root height parameters.
void test_height_transform_jacobian(size_t iter, const char* newick,
                                    int reparameterization, FILE* csv, bool debug) {
    struct timespec start, end;

    Tree* tree = new_Tree(newick, NULL, true);
    parse_dates(tree);
    init_leaf_heights_from_times(tree);

    init_heights_from_bls(tree);

    Model* mtree = new_TreeModel("letree", tree);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_RATIO);
    Model* mtt = mtree->data;
    TreeTransform* tt = mtt->obj;
    Parameters* reparams = get_reparams(tree);
    tt->update(tt);  // update once

    double logP;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        logP = tt->log_jacobian(tt);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);
    if (csv != NULL)
        fprintf(csv, "ratio_transform_jacobian%s,evaluation,off,%f,%f\n",
                (reparameterization == TREE_TRANSFORM_RATIO_NAIVE ? "2" : ""),
                mseconds(start, end) / 1000., logP);

    if (debug) {
        printf("logP %f\n", logP);
    }

    double* gradient = malloc((Tree_tip_count(tree) - 1) * sizeof(double));

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        tt->log_jacobian_gradient(tt, gradient);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));
    if (csv != NULL)
        fprintf(csv, "ratio_transform_jacobian%s,gradient,off,%f,\n",
                (reparameterization == TREE_TRANSFORM_RATIO_NAIVE ? "2" : ""),
                mseconds(start, end) / 1000.);

    if (debug) {
        memset(gradient, 0, (Tree_tip_count(tree) - 1) * sizeof(double));
        tt->log_jacobian_gradient(tt, gradient);
        for (int j = 0; j < Tree_tip_count(tree) - 1; j++) {
            printf("dlogP %f\n", gradient[j]);
        }
    }
    free(gradient);
    mtree->free(mtree);
}

// Transform node ratios to node height.
// Derivatives are wrt ratio/root height parameters.
// At each iteration, node heights are updated from the ratios/root height
// parameters.
void test_height_transform(size_t iter, const char* newick, int reparameterization,
                           FILE* csv, bool debug) {
    struct timespec start, end;

    Tree* tree = new_Tree(newick, NULL, true);
    parse_dates(tree);
    init_leaf_heights_from_times(tree);

    // init_heights_from_distances(tree);
    init_heights_from_bls(tree);

    Model* mtree = new_TreeModel("letree", tree);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_RATIO);
    Model* mtt = mtree->data;
    TreeTransform* tt = mtt->obj;
    Parameters* reparams = get_reparams(tree);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        tt->update(tt);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms\n", iter, mseconds(start, end));
    if (csv != NULL)
        fprintf(csv, "ratio_transform%s,evaluation,off,%f,\n",
                (reparameterization == TREE_TRANSFORM_RATIO_NAIVE ? "2" : ""),
                mseconds(start, end) / 1000.);

    double* gradient = malloc((Tree_tip_count(tree) - 1) * sizeof(double));
    double* height_gradient = malloc((Tree_tip_count(tree) - 1) * sizeof(double));
    for (size_t i = 0; i < Tree_tip_count(tree) - 1; i++) {
        height_gradient[i] = 1.0;
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        tt->update(tt);
        tt->jvp(tt, height_gradient, gradient);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));
    if (csv != NULL)
        fprintf(csv, "ratio_transform%s,gradient,off,%f,\n",
                (reparameterization == TREE_TRANSFORM_RATIO_NAIVE ? "2" : ""),
                mseconds(start, end) / 1000.);

    free(gradient);
    free(height_gradient);
    mtree->free(mtree);
}

// Calculate the constant size coalescent log likelihood.
// The derivatives are wrt population size and the node height parameters
void test_constant(size_t iter, const char* newick, FILE* csv, bool debug) {
    struct timespec start, end;

    Tree* tree = new_Tree(newick, NULL, true);
    parse_dates(tree);

    init_leaf_heights_from_times(tree);

    // init_heights_from_distances(tree);
    init_heights_from_bls(tree);

    Model* mtree = new_TreeModel("letree", tree);

    double value = 4;
    Parameter* N = new_Parameter("theta", value, new_Constraint(0, INFINITY));
    Coalescent* coal = new_ConstantCoalescent(tree, N);

    Model* model = new_CoalescentModel("", coal, mtree);

    double logP;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        coal->need_update_intervals = true;
        logP = model->logP(model);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);
    if (csv != NULL)
        fprintf(csv, "coalescent,evaluation,off,%f,%f\n", mseconds(start, end) / 1000.,
                logP);

    if (debug) {
        printf("logP %f\n", logP);
    }

    Parameters* parameters = new_Parameters(Parameters_count(coal->p));
	Parameters_add_parameters(parameters, coal->p);
	Parameters_add_parameters(parameters, get_reparams(coal->tree));

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        coal->need_update_intervals = true;
        model->gradient(model, parameters);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));
    if (csv != NULL)
        fprintf(csv, "coalescent,gradient,off,%f,\n", mseconds(start, end) / 1000.);

    model->free(model);
    mtree->free(mtree);
    free_Parameters(parameters);
    free_Parameter(N);
}

// Calculate the constant size coalescent log likelihood.
// The derivatives are wrt population size and the node height parameters
void test_skyglide(size_t iter, const char* newick, FILE* csv, bool debug,
                   int intervals, double cutoff) {
    struct timespec start, end;

    Tree* tree = new_Tree(newick, NULL, true);
    parse_dates(tree);

    init_leaf_heights_from_times(tree);

    init_heights_from_distances(tree);

    Model* mtree = new_TreeModel("letree", tree);

    double value = 5.0;
    double* values = dvector(intervals);
    for (size_t i = 0; i < intervals; i++) {
        values[i] = value + i;
    }
    Parameter* popSizes = new_Parameter2("theta", values, intervals, new_Constraint(0, INFINITY));
    free(values);

    Coalescent* coal = new_PiecewiseLinearGridCoalescent(tree, popSizes, intervals,
                                                         cutoff);
    free_Parameter(popSizes);

    Model* model = new_CoalescentModel("", coal, mtree);

    double logP;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        coal->need_update_intervals = true;
        logP = model->logP(model);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);
    if (csv != NULL)
        fprintf(csv, "skyglide,evaluation,off,%f,%f\n", mseconds(start, end) / 1000.,
                logP);

    if (debug) {
        printf("logP %f\n", logP);
    }

    // size_t gradient_size = Coalescent_initialize_gradient(
    //     model, GRADIENT_FLAG_COALESCENT_THETA | GRADIENT_FLAG_TREE_HEIGHTS);
    Parameters* parameters = new_Parameters(10);
    Parameters_add_parameters(parameters, coal->p);
    for(size_t i = 0; i < Tree_node_count(tree); i++){
        Node* node = Tree_node(tree, i);
        if(!Node_isleaf(node)){
            Parameters_add(parameters, node->height);
        }
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        coal->need_update_intervals = true;
        model->gradient(model, parameters);
        // Coalescent_gradient(model);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu skyglide evaluations: %f ms\n", iter, mseconds(start, end));
    if (csv != NULL)
        fprintf(csv, "skyglide,gradient,off,%f,\n", mseconds(start, end) / 1000.);

    model->free(model);
    mtree->free(mtree);
}

// Calculate the log tree likelihood of a time tree (log det Jacbian term is
// included). Derivatives are wrt branch ratios/node height parameters. At each
// iteration, the probability matrices are updated and the whole tree likelihood
// is calculated.
void test_tree_likelihood_time(size_t iter, char* fasta_file, const char* newick,
                               int reparameterization, FILE* csv, bool debug) {
    struct timespec start, end;

    Hashtable* hash2 = new_Hashtable_string(100);
    hashtable_set_key_ownership(hash2, false);
    hashtable_set_value_ownership(hash2, false);

    json_node* json = create_json_tree(json_jc69_strict);

    Sequences* sequences = readSequences(fasta_file);
    sequences->datatype = new_NucleotideDataType();
    SitePattern* patterns = new_SitePattern(sequences);
    free_Sequences(sequences);
    Hashtable_add(hash2, "patterns", patterns);

    Tree* tree = new_Tree(newick, NULL, true);
    parse_dates(tree);
    init_leaf_heights_from_times(tree);

    init_heights_from_distances(tree);

    Model* mtree = new_TreeModel("tree", tree);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_RATIO);
    Hashtable_add(hash2, "tree", mtree);

    Model* mlike = new_TreeLikelihoodModel_from_json(json, hash2);
    SingleTreeLikelihood* tlk = mlike->obj;
    Model* mtt = mtree->data;  // node height transform
    mtree->free(mtree);

    double logP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        mtt->update(mtt, NULL, NULL,
                    0);  // force recalculation of node heights from ratios
        logP = mlike->logP(mlike);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);

    if (debug) {
        printf("logP %f\n", logP);
    }

    Parameters* ps = new_Parameters(1);
    Parameters* reparams = get_reparams(tree);
    Parameters_add_parameters(ps, reparams);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        mtt->update(mtt, NULL, NULL, 0);
        mlike->gradient(mlike, ps);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));

    if (debug) {
        Parameters_zero_grad(ps);
        mlike->gradient(mlike, ps);
        for (int j = 0; j < Parameters_count(ps); j++) {
            Parameter* p = Parameters_at(ps, j);
            for (size_t k = 0; k < Parameter_size(p); k++) {
                printf("dlogP %s %f\n", p->name, p->grad[k]);
            }
        }
    }

    free_Parameters(ps);
    free_SitePattern(patterns);
    mlike->free(mlike);
    free_Hashtable(hash2);
    json_free_tree(json);
}

// Calculate the log tree likelihood of an unrooted tree.
// Derivatives are wrt branch lengths.
// At each iteration, the probability matrices are updated and the whole tree
// likelihood is calculated.
void test_tree_likelihood_unrooted(size_t iter, const char* json_model,
                                   const char* fasta_file, const char* newick,
                                   double scaler, FILE* csv, bool debug) {
    struct timespec start, end;

    Hashtable* hash2 = new_Hashtable_string(100);
    hashtable_set_key_ownership(hash2, false);
    hashtable_set_value_ownership(hash2, false);

    json_node* json = create_json_tree(json_model);

    Sequences* sequences = readSequences(fasta_file);
    sequences->datatype = new_NucleotideDataType();
    SitePattern* patterns = new_SitePattern(sequences);
    free_Sequences(sequences);
    Hashtable_add(hash2, "patterns", patterns);

    Tree* tree = new_Tree(newick, NULL, false);
    Tree_init_branch_lengths(tree);
    Tree_scale_distance(tree, scaler);
    Node* root = Tree_root(tree);
    Parameter* distances = root->distance;
    for (size_t i = 0; i < Parameter_size(distances); i++) {
        if (Parameter_value_at(distances, i) < 1.e-6) {
            Parameter_set_value_at(distances, 1.e-6, i);
        }
    }

    Model* mtree = new_TreeModel("tree", tree);
    Hashtable_add(hash2, "tree", mtree);

    Model* mlike = new_TreeLikelihoodModel_from_json(json, hash2);
    SingleTreeLikelihood* tlk = mlike->obj;
    double logP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        SingleTreeLikelihood_update_all_nodes(tlk);
        tlk->m->need_update = true;
        logP = mlike->logP(mlike);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);
    if (csv != NULL)
        fprintf(csv, "treelikelihood%s,evaluation,off,%f,%f\n",
                (tlk->m->modeltype == GTR ? "GTR" : "JC69"),
                mseconds(start, end) / 1000., logP);
    // logP does not match tree time value because it does not include log det
    // Jacobian

    if (debug) {
        printf("logP %f\n", logP);
    }

    Parameters* ps = new_Parameters(1);
    Parameters_add(ps, root->distance);
    if (tlk->m->modeltype == GTR) {
        Parameters_add(ps, tlk->m->rates_simplex->transform->parameter);
        Parameters_add(ps, tlk->m->simplex->transform->parameter);
    }
    // mlike->prepare_gradient(mlike, ps);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        SingleTreeLikelihood_update_all_nodes(tlk);
        tlk->m->need_update = true;
        mlike->gradient(mlike, ps);
        // TreeLikelihood_gradient(mlike);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));
    if (csv != NULL)
        fprintf(csv, "treelikelihood%s,gradient,off,%f,\n",
                (tlk->m->modeltype == GTR ? "GTR" : "JC69"),
                mseconds(start, end) / 1000.);

    if (debug) {
        Parameters_zero_grad(ps);
        mlike->gradient(mlike, ps);
        for (int j = 0; j < Parameters_count(ps); j++) {
            Parameter* p = Parameters_at(ps, j);
            for (size_t k = 0; k < Parameter_size(p); k++) {
                printf("dlogP %s %f\n", p->name, p->grad[k]);
            }
        }
    }

    free_Parameters(ps);
    free_SitePattern(patterns);
    mtree->free(mtree);
    mlike->free(mlike);
    free_Hashtable(hash2);
    json_free_tree(json);
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        printf(
            "USAGE: benchmarking -r iterations -i alignment-file-name -t "
            "tree-file-name\n\n");
        printf("positional arguments:\n");
        printf("  -r iterations             number of iterations\n");
        printf("  -i alignment-file-name    alignment file\n");
        printf("  -t tree-file-name         tree file in newick format\n\n");
        printf("  --debug                   debug mode (print log likelihood)\n\n");
        exit(0);
    }
    size_t iter = 100;
    char* fasta_file = NULL;
    char* newick = NULL;
    bool debug = false;
    double scaler = 1.0;
    int intervals = -1;
    double cutoff = INFINITY;
    FILE* csv = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0) {
            fasta_file = argv[++i];
        } else if (strcmp(argv[i], "-t") == 0) {
            newick = readTree(argv[++i]);
        } else if (strcmp(argv[i], "-r") == 0) {
            iter = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0) {
            csv = fopen(argv[++i], "w");
            fprintf(csv, "function,mode,JIT,time,logprob\n");
        } else if (strcmp(argv[i], "-s") == 0) {
            scaler = atof(argv[++i]);
        } else if (strcmp(argv[i], "--cutoff") == 0) {
            cutoff = atof(argv[++i]);
        } else if (strcmp(argv[i], "--intervals") == 0) {
            intervals = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--debug") == 0) {
            debug = true;
        }
    }

    printf("Height transform log det Jacobian:\n");
    printf("naive:\n");
    test_height_transform_jacobian(iter, newick, TREE_TRANSFORM_RATIO_NAIVE, csv,
                                   debug);
    printf("efficient:\n");
    test_height_transform_jacobian(iter, newick, TREE_TRANSFORM_RATIO, csv, debug);

    printf("Height transform:\n");
    printf("naive:\n");
    test_height_transform(iter, newick, TREE_TRANSFORM_RATIO_NAIVE, csv, debug);
    printf("efficient:\n");
    test_height_transform(iter, newick, TREE_TRANSFORM_RATIO, csv, debug);

    printf("Constant coalescent:\n");
    test_constant(iter, newick, csv, debug);

    if (intervals > 0 && !isinf(cutoff)) {
        printf("Piecewise linear coalescent:\n");
        test_skyglide(iter, newick, csv, debug, intervals, cutoff);
    }

    printf("Tree likelihood unrooted:\n");
    test_tree_likelihood_unrooted(iter, json_jc69_unrooted, fasta_file, newick, scaler,
                                  csv, debug);
    test_tree_likelihood_unrooted(iter, json_gtr_unrooted, fasta_file, newick, scaler,
                                  csv, debug);

    printf("Tree likelihood time tree:\n");
    test_tree_likelihood_time(iter, fasta_file, newick, 0, csv, debug);
    test_tree_likelihood_time(iter, fasta_file, newick, 1, csv, debug);

    free(newick);
    if (csv != NULL) {
        fclose(csv);
    }
}
