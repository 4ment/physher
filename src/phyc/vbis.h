//
//  vbis.h
//  physher
//
//  Created by Mathieu Fourment on 19/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef vbis_h
#define vbis_h

#include <stdio.h>
#include "vb.h"

typedef struct marginal_vb_t{
	Model* model;
	size_t samples;
	double(*calculate)(struct marginal_vb_t*);
	void(*free)(struct marginal_vb_t*);
}marginal_vb_t;

marginal_vb_t* new_Marginal_VB_from_json(json_node* node, Hashtable* hash);

#endif /* vbis_h */
