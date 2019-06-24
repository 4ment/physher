//
//  sbn.h
//  physher
//
//  Created by Mathieu Fourment on 22/6/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef sbn_h
#define sbn_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"


typedef struct SubSplit{
	unsigned int* Y;
	unsigned int* Z;
}SubSplit;

typedef struct SubSplitPair{
	SubSplit* parent;
	SubSplit* child;
	double count;
}SubSplitPair;

typedef struct SubSplitPairVector{
	SubSplitPair** pairs;
	unsigned int size;
	unsigned int capacity;
}SubSplitPairVector;

typedef struct SBN{
	unsigned int** bitsets;
	unsigned int size;
	SubSplit** rootProb;
	double* rootProbWeights;
	unsigned int rootProbCount;
	unsigned int capacity;
	Hashtable* conditionalCladeMapY;
	Hashtable* conditionalCladeMapZ;
}SBN;

SBN* new_SBN_from_json(json_node* node, Hashtable* hash);

#endif /* sbn_h */
