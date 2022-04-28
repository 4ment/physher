#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "phyc/datatype.h"
#include "phyc/demographicmodels.h"
#include "phyc/parameters.h"
#include "phyc/sequenceio.h"
#include "phyc/sitepattern.h"
#include "phyc/tree.h"
#include "phyc/treeio.h"
#include "phyc/treelikelihood.h"
#include "phyc/treetransform.h"

double mseconds(struct timespec start, struct timespec end) {
    return (end.tv_sec - start.tv_sec) * 1000. + (end.tv_nsec - start.tv_nsec) / 1000000.;
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
			\"type\":\"Simplex\", \
			\"values\":[0.25,0.25,0.25,0.25] \
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
			\"id\":\"rate\", \"type\":\"parameter\", \"value\":0.001, \"lower\":0 \
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
			\"type\":\"Simplex\", \
			\"values\":[0.25,0.25,0.25,0.25] \
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
		\"rates\":{ \
			\"ac\":{\"id\":\"ac\", \"type\":\"parameter\", \"value\":1, \"lower\":0, \"upper\":\"infinity\"}, \
			\"ag\":{\"id\":\"ag\", \"type\":\"parameter\", \"value\":1, \"lower\":0, \"upper\":\"infinity\"}, \
			\"at\":{\"id\":\"at\", \"type\":\"parameter\", \"value\":1, \"lower\":0, \"upper\":\"infinity\"}, \
			\"cg\":{\"id\":\"cg\", \"type\":\"parameter\", \"value\":1, \"lower\":0, \"upper\":\"infinity\"}, \
			\"ct\":{\"id\":\"ct\", \"type\":\"parameter\", \"value\":1, \"lower\":0, \"upper\":\"infinity\"} \
		}, \
		\"frequencies\":{ \
			\"id\":\"freqs\", \
			\"type\":\"Simplex\", \
			\"values\":[0.25,0.25,0.25,0.25] \
		} \
	}, \
	\"sitemodel\":{ \
		\"id\": \"sitemodel\", \
		\"type\": \"sitemodel\" \
	 \
	}, \
	\"tree\": \"&tree\" \
}";

// Calculate the log of the determinant of the Jacobian for the node height transform.
// Derivatives are wrt ratio/root height parameters.
void test_height_transform_jacobian(size_t iter, const char* newick, int reparameterization, FILE* csv, bool debug) {
    struct timespec start, end;

    Tree* tree = new_Tree(newick, true);
    Tree_set_transform(tree, reparameterization);
    parse_dates(tree);
    init_leaf_heights_from_times(tree);

    init_heights_from_distances(tree);

    Model* mtree = new_TreeModel("letree", tree);
    Model* mtt = mtree->data;
    TreeTransform* tt = mtt->obj;
    Parameters* reparams = get_reparams(tree);
    tt->update(tt);  // update once

    double logP;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        logP = mtt->logP(mtt);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);
    if (csv != NULL)
        fprintf(csv, "ratio_transform_jacobian%s,evaluation,off,%f,%f\n", (reparameterization == 1 ? "2" : ""), mseconds(start, end) / 1000., logP);

    if (debug) {
        printf("logP %f\n", logP);
    }

    double dlogP;
    double* gradient = malloc((Tree_tip_count(tree) - 1) * sizeof(double));

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        tt->log_jacobian_gradient(tt, gradient);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));
    if (csv != NULL)
        fprintf(csv, "ratio_transform_jacobian%s,gradient,off,%f,\n", (reparameterization == 1 ? "2" : ""), mseconds(start, end) / 1000.);

    if (debug) {
        for (int j = 0; j < Parameters_count(reparams); j++) {
            logP = mtt->dlogP(mtt, Parameters_at(reparams, j));
            printf("dlogP %f\n", dlogP);
        }
    }
    free(gradient);
    mtree->free(mtree);
}

// Transform node ratios to node height.
// Derivatives are wrt ratio/root height parameters.
// At each iteration, node heights are updated from the ratios/root height parameters.
void test_height_transform(size_t iter, const char* newick, int reparameterization, FILE* csv, bool debug) {
    struct timespec start, end;

    Tree* tree = new_Tree(newick, true);
    Tree_set_transform(tree, reparameterization);
    parse_dates(tree);
    init_leaf_heights_from_times(tree);

    init_heights_from_distances(tree);

    Model* mtree = new_TreeModel("letree", tree);
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
        fprintf(csv, "ratio_transform%s,evaluation,off,%f,\n", (reparameterization == 1 ? "2" : ""), mseconds(start, end) / 1000.);

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
        fprintf(csv, "ratio_transform%s,gradient,off,%f,\n", (reparameterization == 1 ? "2" : ""), mseconds(start, end) / 1000.);

    free(gradient);
    free(height_gradient);
    mtree->free(mtree);
}

// Calculate the constant size coalescent log likelihood.
// The derivatives are wrt population size and the node height parameters
void test_constant(size_t iter, const char* newick, FILE* csv, bool debug) {
    struct timespec start, end;

    Tree* tree = new_Tree(newick, true);
    parse_dates(tree);

    init_leaf_heights_from_times(tree);

    init_heights_from_distances(tree);

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
        fprintf(csv, "coalescent,evaluation,off,%f,%f\n", mseconds(start, end) / 1000., logP);

    if (debug) {
        printf("logP %f\n", logP);
    }

    Parameters* ps = new_Parameters(1);
    Parameters_add(ps, N);
    Parameter_set_value(N, value);
    Node** nodes = Tree_nodes(tree);
    for (int i = 0; i < Tree_node_count(tree); i++) {
        if (!Node_isleaf(nodes[i])) {
            Parameters_add(ps, nodes[i]->height);
        }
    }

    model->prepare_gradient(model, ps);
    double dlogP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        coal->need_update_intervals = true;
        for (int j = 0; j < Parameters_count(ps); j++) {
            dlogP = model->dlogP(model, Parameters_at(ps, j));
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));
    if (csv != NULL)
        fprintf(csv, "coalescent,gradient,off,%f,\n", mseconds(start, end) / 1000.);

    model->free(model);
    mtree->free(mtree);
    free_Parameter(N);
    free_Parameters(ps);
}

// Calculate the log tree likelihood of a time tree (log det Jacbian term is included).
// Derivatives are wrt branch ratios/node height parameters.
// At each iteration, the probability matrices are updated and the whole tree likelihood
// is calculated.
void test_tree_likelihood_time(size_t iter, char* fasta_file, const char* newick, int reparameterization, FILE* csv, bool debug) {
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

    Tree* tree = new_Tree(newick, true);
    Tree_set_transform(tree, reparameterization);
    parse_dates(tree);
    init_leaf_heights_from_times(tree);

    init_heights_from_distances(tree);

    Model* mtree = new_TreeModel("tree", tree);
    Hashtable_add(hash2, "tree", mtree);

    Model* mlike = new_TreeLikelihoodModel_from_json(json, hash2);
    SingleTreeLikelihood* tlk = mlike->obj;
    Model* mtt = mtree->data;  // node height transform
    mtree->free(mtree);

    double logP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        mtt->update(mtt, NULL, 0);  // force recalculation of node heights from ratios
        logP = mlike->logP(mlike);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);

    if (debug) {
        printf("logP %f\n", logP);
    }

    Parameters* ps = new_Parameters(1);
    Parameters* reparams = get_reparams(tree);
    for (int i = 0; i < Parameters_count(reparams); i++) {
        Parameters_add(ps, Parameters_at(reparams, i));
    }
    mlike->prepare_gradient(mlike, ps);
    double dlogP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        mtt->update(mtt, NULL, 0);
        for (int j = 0; j < Parameters_count(ps); j++) {
            dlogP = mlike->dlogP(mlike, Parameters_at(ps, j));
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));

    if (debug) {
        for (int j = 0; j < Parameters_count(ps); j++) {
            dlogP = mtt->dlogP(mtt, Parameters_at(ps, j));
            printf("dlogP %f\n", dlogP);
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
// At each iteration, the probability matrices are updated and the whole tree likelihood
// is calculated.
void test_tree_likelihood_unrooted(size_t iter, const char* json_model, const char* fasta_file, const char* newick, double scaler, FILE* csv, bool debug) {
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

    Tree* tree = new_Tree(newick, false);
    Tree_scale_distance(tree, scaler);
    Node* root = Tree_root(tree);
    for (int i = 0; i < Tree_node_count(tree); i++) {
        Node* node = Tree_node(tree, i);
        if (node != root) {
            if (Parameter_value(node->distance) < 1.e-6) {
                Parameter_set_value(node->distance, 1.e-6);
            }
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
        fprintf(csv, "treelikelihood%s,evaluation,off,%f,%f\n", (tlk->m->modeltype == GTR ? "GTR" : "JC69"), mseconds(start, end) / 1000., logP);
    // logP does not match tree time value because it does not include log det Jacobian

    if (debug) {
        printf("logP %f\n", logP);
    }

    Parameters* ps = new_Parameters(1);
    for (int i = 0; i < Tree_node_count(tree); i++) {
        Node* node = Tree_node(tree, i);
        if (node != root && root->right != node) {
            Parameters_add(ps, node->distance);
        }
    }
    if (tlk->m->modeltype == GTR) {
        Parameters_add_parameters(ps, tlk->m->rates);
        Parameters_add_parameters(ps, tlk->m->simplex->parameters);
    }
    mlike->prepare_gradient(mlike, ps);
    double dlogP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        SingleTreeLikelihood_update_all_nodes(tlk);
        tlk->m->need_update = true;
        for (int j = 0; j < Parameters_count(ps); j++) {
            dlogP = mlike->dlogP(mlike, Parameters_at(ps, j));
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));
    if (csv != NULL)
        fprintf(csv, "treelikelihood%s,gradient,off,%f,\n", (tlk->m->modeltype == GTR ? "GTR" : "JC69"), mseconds(start, end) / 1000.);

    if (debug) {
        for (int j = 0; j < Parameters_count(ps); j++) {
            dlogP = mlike->dlogP(mlike, Parameters_at(ps, j));
            printf("dlogP %f\n", dlogP);
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
        printf("USAGE: benchmarking -r iterations -i alignment-file-name -t tree-file-name\n\n");
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
        } else if (strcmp(argv[i], "--debug") == 0) {
            debug = true;
        }
    }

    printf("Height transform log det Jacobian:\n");
    printf("naive:\n");
    test_height_transform_jacobian(iter, newick, 0, csv, debug);
    printf("efficient:\n");
    test_height_transform_jacobian(iter, newick, 1, csv, debug);

    printf("Height transform:\n");
    printf("naive:\n");
    test_height_transform(iter, newick, 0, csv, debug);
    printf("efficient:\n");
    test_height_transform(iter, newick, 1, csv, debug);

    printf("Constant coalescent:\n");
    test_constant(iter, newick, csv, debug);

    printf("Tree likelihood unrooted:\n");
    test_tree_likelihood_unrooted(iter, json_jc69_unrooted, fasta_file, newick, scaler, csv, debug);
    test_tree_likelihood_unrooted(iter, json_gtr_unrooted, fasta_file, newick, scaler, csv, debug);

    printf("Tree likelihood time tree:\n");
    test_tree_likelihood_time(iter, fasta_file, newick, 0, csv, debug);
    test_tree_likelihood_time(iter, fasta_file, newick, 1, csv, debug);

    free(newick);
    if (csv != NULL) {
        fclose(csv);
    }
}
