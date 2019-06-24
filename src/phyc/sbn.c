//
//  sbn.c
//  physher
//
//  Created by Mathieu Fourment on 22/6/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "sbn.h"

#include "treeio.h"
#include "hashfunctions.h"
#include "matrix.h"

#define WORD_BITS (8 * sizeof(unsigned int))

#define PRECOMPUTE 1

unsigned int HashBitset(const void* data) {
	const unsigned int* bitset = data;
	unsigned int hash = bitset[0];
	for(int i = 1; i <= bitset[0]; i++){
		hash ^= bitset[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
	}
	return hash;
}

unsigned int preComputedHashBitset(const void* data) {
	const unsigned int* bitset = data;
	return bitset[1];
}

void setHashBitset(unsigned int* bitset) {
	unsigned int hash = bitset[0];
	for(int i = 2; i <= bitset[0]; i++){
		hash ^= bitset[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
	}
	bitset[1] = hash;
}

// From java Arrays.hasCode
unsigned int HashSubSplit(const void* data) {
	const SubSplit* subSplit = data;
#ifdef PRECOMPUTE
	return 31*(31 + subSplit->Y[1]) + subSplit->Z[1];
#else
	return 31*(31 + HashBitset(subSplit->Y)) + HashBitset(subSplit->Z);
#endif
}

SubSplitPair* new_SubSplitPair(SubSplit* parent, SubSplit* child){
	SubSplitPair* pair = malloc(sizeof(SubSplitPair));
	pair->parent = parent;
	pair->child = child;
	pair->count = 0;
	return pair;
}

bool SubSplit_equals(const SubSplit* a, const SubSplit* b){
#ifdef PRECOMPUTE
	if(a->Y[1] != b->Y[1] || a->Z[1] != b->Z[1]) return false;
	if(memcmp(a->Y+2, b->Y+2, a->Y[0]*sizeof(unsigned int)) != 0) return false;
	return memcmp(a->Z+2, b->Z+2, a->Y[0]*sizeof(unsigned int)) == 0;
//	unsigned int count = a->Y[0]+1;
//	for (int i = 2; i <= count; i++) {
//		if (a->Y[i] != b->Y[i] || a->Z[i] != b->Z[i]) return false;
//	}
#else
	if(memcmp(a->Y+1, b->Y+1, a->Y[0]*sizeof(unsigned int)) != 0) return false;
	return memcmp(a->Z+1, b->Z+1, a->Y[0]*sizeof(unsigned int)) == 0;
#endif
	return true;
}

bool SubSplitPair_equals(SubSplitPair* a, SubSplitPair* b){
	if(!SubSplit_equals(a->parent, b->parent)) return false;
	return SubSplit_equals(a->child, b->child);
}

static bool Hashtable_compare_SubSplit( const void *key1, const void *key2 ){
	return SubSplit_equals(key1, key2);
}


static bool Hashtable_compare_bitset( const void *key1, const void *key2 ){
	const unsigned int* k1 = (unsigned int*)key1;
	const unsigned int* k2 = (unsigned int*)key2;
#ifdef PRECOMPUTE
	return memcmp(k1+2, k2+2, k1[0]*sizeof(unsigned int)) == 0;
#else
	return memcmp(k1+1, k2+1, k1[0]*sizeof(unsigned int)) == 0;
#endif
}


void addPair(Hashtable* conditionalCladeMap, SubSplitPair* pair){
	void* pairs = Hashtable_get(conditionalCladeMap, pair->parent);
	
	// subSplitPair does not exists
	if (pairs == NULL){
		SubSplitPairVector* vec = malloc(sizeof(SubSplitPairVector));
		vec->capacity = 10;
		vec->pairs = malloc(vec->capacity*sizeof(SubSplitPair*));
		vec->size = 1;
		vec->pairs[0] = pair;
		pair->count++;
		Hashtable_add(conditionalCladeMap, pair->parent, vec);
	}
	else {
		SubSplitPairVector* vec = pairs;
		int i = 0;
		for ( ; i < vec->size; i++){
			// found and increment
			if (SubSplitPair_equals(vec->pairs[i], pair)){
				vec->pairs[i]->count++;
				free(pair);
				break;
			}
		}
		// parent subsplit exists but not the child subsplit
		if (vec->size == i){
			if(vec->size == vec->capacity){
				vec->capacity += 10;
				vec->pairs = realloc(vec->pairs, vec->capacity*sizeof(SubSplitPair*));
			}
			vec->pairs[vec->size] = pair;
			pair->count++;
			vec->size++;
		}
	}
}

static inline void bitset_or(unsigned int* bitarray, unsigned int* bitarray1, unsigned int* bitarray2) {
	bitarray[0] = bitarray1[0];
#ifdef PRECOMPUTE
	size_t len = bitarray[0]+1;
	for (int i = 2; i <= len; i++) {
#else
	for (int i = 1; i <= bitarray[0]; i++) {
#endif
		bitarray[i] = bitarray1[i] | bitarray2[i];
	}
}

static inline void bitset_set(unsigned int* bitarray, size_t idx) {
#ifdef PRECOMPUTE
	bitarray[2+idx/WORD_BITS] |= (1 << (idx % WORD_BITS));
#else
	bitarray[1+idx/WORD_BITS] |= (1 << (idx % WORD_BITS));
#endif
}

static inline unsigned int bitset_next_set_bit(unsigned int* bitset){
	size_t len = bitset[0];
#ifdef PRECOMPUTE
	len++;
	for (int i = 2; i <= len; i++) {
		if (bitset[i] != 0) {
			for (int j = 0; j < WORD_BITS; j++) {
				if((bitset[i] >> j) & 1U) return (i-2)*WORD_BITS+j;
			}
		}
	}
#else
	for (int i = 1; i <= len; i++) {
		if (bitset[i] != 0) {
			for (int j = 0; j < WORD_BITS; j++) {
				if((bitset[i] >> j) & 1U) return (i-1)*WORD_BITS+j;
			}
		}
	}
#endif
	return -1;
}

static inline bool bitset_equals(unsigned int* bitset1, unsigned int* bitset2){
#ifdef PRECOMPUTE
	//	if(bitset1[1] != bitset2[1]) return false;
	return memcmp(bitset1+2, bitset2+2, bitset1[0]*sizeof(unsigned int)) == 0;
#else
	return memcmp(bitset1+1, bitset2+1, bitset1[0]*sizeof(unsigned int)) == 0;
#endif
}

void process_tree(SBN* sbn, Node* node, SubSplit* ssParent){
	if (!Node_isleaf(node)) {
		unsigned int* Y = sbn->bitsets[Node_id(node->left)];
		unsigned int* Z = sbn->bitsets[Node_id(node->right)];
		if (bitset_next_set_bit(Y) < bitset_next_set_bit(Z)){
			unsigned int* temp = Y;
			Y = Z;
			Z = temp;
		}
		SubSplit* ss = malloc(sizeof(SubSplit));
		ss->Y = Y;
		ss->Z = Z;
		
		if (Node_isroot(node)) {
			int index = 0;
			for(; index < sbn->rootProbCount; index++){
				if (SubSplit_equals(sbn->rootProb[index], ss)) {
					free(ss);
					ss = sbn->rootProb[index];
					break;
				}
			}
			if (index == sbn->rootProbCount){
				if(sbn->rootProbCount == sbn->capacity && sbn->rootProbCount > 0){
					sbn->capacity += 10;
					sbn->rootProb = realloc(sbn->rootProb, sbn->capacity*sizeof(SubSplit*));
					sbn->rootProbWeights = realloc(sbn->rootProbWeights, sbn->capacity*sizeof(unsigned int));
					memset(sbn->rootProbWeights+sbn->rootProbCount, 0, 10*sizeof(unsigned int));
				}
				sbn->rootProb[sbn->rootProbCount] = ss;
				index = sbn->rootProbCount++;
			}
			sbn->rootProbWeights[index]++;
		} else {
			unsigned int* X = sbn->bitsets[Node_id(node)];
			SubSplitPair* pair = new_SubSplitPair(ssParent, ss);
			if (bitset_equals(X, ssParent->Y)) {
				addPair(sbn->conditionalCladeMapY, pair);
			}
			else{
				addPair(sbn->conditionalCladeMapZ, pair);
			}
		}
		process_tree(sbn, node->left, ss);
		process_tree(sbn, node->right, ss);
	}
}

void normalizeCPT_aux(Hashtable* conditionalCladeMap){
	double total;
	Hashtable_init_iterator(conditionalCladeMap);
	HashEntry *entry = NULL;
	while ( (entry = Hashtable_next(conditionalCladeMap) ) != NULL ) {
		SubSplitPairVector *vec = (SubSplitPairVector*)HashEntry_value(entry);
		total = 0;
		for (int i = 0; i < vec->size; i++) {
			total += vec->pairs[i]->count;
		}
		for (int i = 0; i < vec->size; i++) {
			vec->pairs[i]->count /= total;
		}
	}
}

void normalizeCPT(SBN* sbn){
	double total = 0;
	for (int i = 0; i < sbn->rootProbCount; i++){
		total += sbn->rootProbWeights[i];
	}
	for (int i = 0; i < sbn->rootProbCount; i++){
		sbn->rootProbWeights[i] /= total;
	}

	normalizeCPT_aux(sbn->conditionalCladeMapY);
	normalizeCPT_aux(sbn->conditionalCladeMapZ);
}
	
SBN* build(const char* file){
	SBN* sbn = malloc(sizeof(SBN));
	sbn->bitsets = NULL;
	sbn->conditionalCladeMapY = new_Hashtable( 100, &HashSubSplit, Hashtable_compare_SubSplit, HASHTABLE_KEY_REFERENCE);
	sbn->conditionalCladeMapZ = new_Hashtable( 100, &HashSubSplit, Hashtable_compare_SubSplit, HASHTABLE_KEY_REFERENCE);
	hashtable_set_key_ownership(sbn->conditionalCladeMapY, false);
	hashtable_set_value_ownership(sbn->conditionalCladeMapY, false);
	hashtable_set_key_ownership(sbn->conditionalCladeMapZ, false);
	hashtable_set_value_ownership(sbn->conditionalCladeMapZ, false);
	sbn->capacity = 10;
	sbn->rootProbCount = 0;
	sbn->rootProb = malloc(sbn->capacity*sizeof(SubSplit*));
	sbn->rootProbWeights = calloc(sbn->capacity,sizeof(double));
	Hashtable* hash = NULL;
#ifdef PRECOMPUTE
	Hashtable* bitsets = new_Hashtable( 100, preComputedHashBitset, Hashtable_compare_bitset, HASHTABLE_KEY_REFERENCE);
#else
	Hashtable* bitsets = new_Hashtable( 100, HashBitset, Hashtable_compare_bitset, HASHTABLE_KEY_REFERENCE);
#endif
	hashtable_set_key_ownership(bitsets, false);
	hashtable_set_value_ownership(bitsets, false);
	int* indexes = NULL;
	TreeFileIterator *iter = new_TreeFileIterator(file);
	char *line = NULL;
	unsigned int* temp = NULL;
	unsigned int bitset_length;

	while ( (line = TreeFileIterator_next_tree(iter)) != NULL  ) {
		Tree *t = new_Tree(line, true);
		free(line);
		Node** nodes = Tree_get_nodes(t, POSTORDER);
		if (sbn->bitsets == NULL) {
			sbn->size = Tree_node_count(t);
			sbn->bitsets = malloc(sbn->size*sizeof(unsigned int*));
#ifdef PRECOMPUTE
			bitset_length = (sbn->size+WORD_BITS-1)/WORD_BITS + 2;
			temp = calloc(bitset_length, sizeof(unsigned int));
			temp[0] = bitset_length-2;
#else
			bitset_length = (sbn->size+WORD_BITS-1)/WORD_BITS + 1;
			temp = calloc(bitset_length, sizeof(unsigned int));
			temp[0] = bitset_length-1;
#endif
			indexes = ivector(Tree_tip_count(t));
			hash = new_Hashtable_string(Tree_tip_count(t));
			hashtable_set_key_ownership(hash, true);
			hashtable_set_value_ownership(hash, false);
			for (int i = 0; i < Tree_tip_count(t); i++) {
				indexes[i] = i;
			}
			for (int i = 0; i < sbn->size; i++) {
				if(Node_isleaf(nodes[i])){
					Hashtable_add(hash, String_clone(Node_name(nodes[i])), indexes+Node_class_id(nodes[i]));
				}
			}
		}

		int internal = Tree_tip_count(t);
		for (int i = 0; i < sbn->size; i++) {
			if(Node_isleaf(nodes[i])){
				int index = *((int*)Hashtable_get(hash, Node_name(nodes[i])));
				Node_set_id(nodes[i], index);
			}
			else{
				Node_set_id(nodes[i], internal++);
			}
		}

		for (int i = 0; i < sbn->size; i++) {
			int index = Node_id(nodes[i]);
			if (Node_isleaf(nodes[i])) {
				memset(temp+1, 0, (bitset_length-1)*sizeof(unsigned int));
				bitset_set(temp, index);
			}
			else{
				bitset_or(temp, sbn->bitsets[Node_id(Node_left(nodes[i]))], sbn->bitsets[Node_id(Node_right(nodes[i]))]);
			}
#ifdef PRECOMPUTE
			setHashBitset(temp);
#endif

			void* data;
			if ((data = Hashtable_get(bitsets, temp)) != NULL) {
				sbn->bitsets[index] = data;
			}
			else{
				unsigned int* bitset = calloc(bitset_length, sizeof(unsigned int));
				memcpy(bitset, temp, bitset_length*sizeof(unsigned int));
				Hashtable_add(bitsets, bitset, bitset);
				sbn->bitsets[index] = bitset;
			}
		}
		process_tree(sbn, Tree_root(t), NULL);
		free_Tree(t);
	}
	printf("%d\n", Hashtable_length(sbn->conditionalCladeMapY));
	printf("%d\n", Hashtable_length(sbn->conditionalCladeMapZ));
	free(temp);
	free(indexes);
	free_Hashtable(hash);
	free_Hashtable(bitsets);
	free_TreeFileIterator(iter);
	
	normalizeCPT(sbn);
	for (int i = 0; i < sbn->rootProbCount; i++){
		printf("%f\n", sbn->rootProbWeights[i]);
	}
	return sbn;
}

void free_SBN(SBN* sbn){
	for (int i = 0; i < sbn->size; i++) free(sbn->bitsets[i]);
	free(sbn->bitsets);
	free(sbn->rootProbWeights);
	free_Hashtable(sbn->conditionalCladeMapY);
	free_Hashtable(sbn->conditionalCladeMapZ);
	free(sbn);
}
	
SBN* new_SBN_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"trees"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));

	char* filepath = get_json_node_value_string(node, "trees");
	return build(filepath);
}
