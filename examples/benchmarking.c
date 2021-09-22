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
	\"sitemodel\":{ \
		\"id\": \"sitemodel\", \
		\"type\": \"sitemodel\", \
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
		} \
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
	\"sitemodel\":{ \
		\"id\": \"sitemodel\", \
		\"type\": \"sitemodel\", \
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
		} \
	 \
	}, \
	\"tree\": \"&tree\" \
}";

// Calculate the log of the determinant of the Jacobian for the node height transform.
// Derivatives are wrt ratio/root height parameters.
// At each iteration, node heights are updated from the ratios/root height parameters.
void test_height_transform(size_t iter, const char* newick, bool debug) {
    struct timespec start, end;

    Tree* tree = new_Tree(newick, true);
    parse_dates(tree);
    init_leaf_heights_from_times(tree);

    init_heights_from_distances(tree);

    Model* mtree = new_TreeModel("letree", tree);
    Model* mtt = mtree->data;
    TreeTransform* tt = mtt->obj;
    Parameters* reparams = get_reparams(tree);

    double logP;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        tt->update(tt);
        logP = mtt->logP(mtt);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);

    if (debug) {
        printf("logP %f\n", logP);
    }

    double dlogP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        tt->update(tt);
        for (int j = 0; j < Parameters_count(reparams); j++) {
            logP = mtt->dlogP(mtt, Parameters_at(reparams, j));
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));

    if (debug) {
        for (int j = 0; j < Parameters_count(reparams); j++) {
            logP = mtt->dlogP(mtt, Parameters_at(reparams, j));
            printf("dlogP %f\n", dlogP);
        }
    }

    mtree->free(mtree);
}

// Calculate the constant size coalescent log likelihood.
// The derivatives are wrt population size and the ratio/node height parameters
// but the log determinant of the Jacobian is not included.
// At each iteration, node heights are updated from the ratios/root height parameters.
void test_constant(size_t iter, const char* newick, bool debug) {
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

    Model* mtt = mtree->data;
    double logP;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        //      	coal->need_update_intervals = true;
        mtt->update(mtt, NULL, 0);
        logP = model->logP(model);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);

    if (debug) {
        printf("logP %f\n", logP);
    }

    Parameters* ps = new_Parameters(1);
    Parameters_add(ps, N);
    Parameter_set_value(N, value);
    Parameters* reparams = get_reparams(tree);
    for (int i = 0; i < Parameters_count(reparams); i++) {
        Parameters_add(ps, Parameters_at(reparams, i));
    }

    model->prepare_gradient(model, ps);
    double dlogP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        coal->need_update_intervals = true;
        mtt->update(mtt, NULL, 0);
        for (int j = 0; j < Parameters_count(ps); j++) {
            dlogP = model->dlogP(model, Parameters_at(ps, j));
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

    model->free(model);
    mtree->free(mtree);
    free_Parameter(N);
    free_Parameters(ps);
}

// Calculate the log tree likelihood of a time tree (log det Jacbian term is included).
// Derivatives are wrt branch ratios/node height parameters.
// At each iteration, the probability matrices are updated and the whole tree likelihood
// is calculated.
void test_tree_likelihood_time(size_t iter, char* fasta_file, const char* newick, bool debug) {
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
    parse_dates(tree);
    init_leaf_heights_from_times(tree);

    init_heights_from_distances(tree);

    Model* mtree = new_TreeModel("tree", tree);
    Hashtable_add(hash2, "tree", mtree);

    Model* mlike = new_TreeLikelihoodModel_from_json(json, hash2);
    SingleTreeLikelihood* tlk = mlike->obj;
    Model* mtt = mtree->data;  // node height transform

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
void test_tree_likelihood_unrooted(size_t iter, const char* fasta_file, const char* newick, bool debug) {
    struct timespec start, end;

    Hashtable* hash2 = new_Hashtable_string(100);
    hashtable_set_key_ownership(hash2, false);
    hashtable_set_value_ownership(hash2, false);

    json_node* json = create_json_tree(json_jc69_unrooted);

    Sequences* sequences = readSequences(fasta_file);
    sequences->datatype = new_NucleotideDataType();
    SitePattern* patterns = new_SitePattern(sequences);
    free_Sequences(sequences);
    Hashtable_add(hash2, "patterns", patterns);

    Tree* tree = new_Tree(newick, false);
    Tree_scale_distance(tree, 0.001);
    Model* mtree = new_TreeModel("tree", tree);
    Hashtable_add(hash2, "tree", mtree);

    Model* mlike = new_TreeLikelihoodModel_from_json(json, hash2);
    SingleTreeLikelihood* tlk = mlike->obj;
    double logP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        SingleTreeLikelihood_update_all_nodes(tlk);
        logP = mlike->logP(mlike);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu evaluations: %f ms (%f)\n", iter, mseconds(start, end), logP);
    // logP does not match tree time value because it does not include log det Jacobian

    if (debug) {
        printf("logP %f\n", logP);
    }

    Parameters* ps = new_Parameters(1);
    Node* root = Tree_root(tree);
    for (int i = 0; i < Tree_node_count(tree); i++) {
        Node* node = Tree_node(tree, i);
        if (node != root && root->right != node) {
            Parameters_add(ps, node->distance);
        }
    }
    mlike->prepare_gradient(mlike, ps);
    double dlogP;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (size_t i = 0; i < iter; i++) {
        SingleTreeLikelihood_update_all_nodes(tlk);
        for (int j = 0; j < Parameters_count(ps); j++) {
            dlogP = mlike->dlogP(mlike, Parameters_at(ps, j));
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("  %zu gradient evaluations: %f ms\n", iter, mseconds(start, end));

    if (debug) {
        for (int j = 0; j < Parameters_count(ps); j++) {
            dlogP = mlike->dlogP(mlike, Parameters_at(ps, j));
            printf("dlogP %f\n", dlogP);
        }
    }

    free_Parameters(ps);
    free_SitePattern(patterns);
    mlike->free(mlike);
    free_Hashtable(hash2);
    json_free_tree(json);
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        printf("USAGE: benchmarking iterations alignment-file-name tree-file-name\n\n");
        printf("positional arguments:\n");
        printf("  iterations             number of iterations\n");
        printf("  alignment-file-name    alignment file\n");
        printf("  tree-file-name number  tree file in newick format\n\n");
        exit(0);
    }
    size_t iter = atoi(argv[1]);
    char* fasta_file = argv[2];
    char* newick = readTree(argv[3]);
    bool debug = argc > 4 ? atoi(argv[4]) : false;

    printf("Height transform log det Jacobian:\n");
    test_height_transform(iter, newick, debug);

    printf("Constant coalescent:\n");
    test_constant(iter, newick, debug);

    printf("Tree likelihood unrooted:\n");
    test_tree_likelihood_unrooted(iter, fasta_file, newick, debug);

    printf("Tree likelihood time tree:\n");
    test_tree_likelihood_time(iter, fasta_file, newick, debug);
}