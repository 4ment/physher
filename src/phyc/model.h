/*
 *  model.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/15/10.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef _MODEL_H_
#define _MODEL_H_

#include <stdio.h>

#include "utils.h"
#include "mjson.h"
#include "hashtable.h"

struct _Parameter;
typedef struct _Parameter Parameter;

struct _Parameters;
typedef struct _Parameters Parameters;

struct _Constraint;
typedef struct _Constraint Constraint;

struct _ListenerList;
typedef struct _ListenerList ListenerList;

struct _Listeners;
typedef struct _Listener Listener;

struct _Model;
typedef struct _Model Model;

typedef enum model_t{
	MODEL_ALIGNMENT=0,
	MODEL_BOUNDMODEL,
	MODEL_BRANCHMODEL,
	MODEL_COALESCENT,
	MODEL_COMPOUND,
	MODEL_DISCRETE_PARAMETER,
	MODEL_DISTRIBUTION,
	MODEL_JACOBIAN_TRANSFORM,
	MODEL_LAPLACE,
	MODEL_PARAMETERS,
	MODEL_PARSIMONY,
	MODEL_SITEMODEL,
	MODEL_SUBSTITUTION,
	MODEL_TREE,
	MODEL_TREE_TRANSFORM,
	MODEL_TREELIKELIHOOD,
    MODEL_VARIATIONAL
}model_t;

static const char* model_type_strings[] = {
	"alignment",
	"bound",
	"branchmodel",
	"coalescent",
	"compound",
	"discreteparameter",
	"distribution",
	"jacobiantransform",
	"laplace",
	"parameters",
	"parsimony",
	"sitemodel",
	"substitutionmodel",
	"tree",
	"treetransform",
	"treelikelihood",
    "variational"
};

model_t check_model(const char* type);

#pragma mark -
#pragma mark ListenerList

struct _ListenerList {
    Model **models;
    int count;
    int capacity;
    bool enabled;
    Parameters *parameters;
    void (*free)(ListenerList *);
    void (*fire)(ListenerList *, Model *, Parameter *, int);
    void (*add)(ListenerList *, Model *);
	void (*add_parameter)(ListenerList *, Parameter *);
    void (*remove)(ListenerList *, Model *);
    void (*removeAll)(ListenerList *);
};

ListenerList * new_ListenerList( const unsigned capacity );

#pragma mark -
#pragma mark Model

// Layout requested from a Model's hessian function. The matrix is taken over
// the flattened scalar elements of the Parameters argument (each Parameter
// expanded by Parameter_size, in list order). dim = that total element count.
//  - HESSIAN_DIAGONAL: out has dim entries,   out[k]        = d2 logP / dx_k^2
//  - HESSIAN_FULL:     out has dim*dim entries (row-major, symmetric),
//                      out[i*dim + j] = d2 logP / dx_i dx_j
// Derivatives are in the natural (constrained) parameter space; any transform
// to unconstrained space is the caller's responsibility.
typedef enum {
	HESSIAN_DIAGONAL = 0,
	HESSIAN_FULL = 1,
} hessian_mode_t;

struct _Model {
	void *obj; // pointer to model
	char *name;
	model_t type;
	void* data;
	double (*logP)( Model * );
	double (*full_logP)( Model * );
	void (*gradient)( Model *, Parameters*);
	void (*hessian)( Model *, const Parameters*, hessian_mode_t, double*);
	Model* (*clone)( Model *, Hashtable* );
	void (*free)( Model * );
	void (*update)( Model *, Model *, Parameter*, int );
	void (*reset)(Model*);
	void (*sample)(Model *);
	void (*rsample)(Model *);

	// if the model is a transform like a simplex this function retrieve/set the
	// values
	void (*get)(Model *, double *);
	void (*set)(Model *, const double *);

    ListenerList *listeners;
	int ref_count;
    Parameters *parameters;  // parameters of the model, including from submodels

    void(*store)(Model*);
	void(*restore)(Model*);
	void(*accept)(Model*);
	double lp;
	double storedLogP;
	double stored;
	bool samplable; // model is a distribution that can sampled directly
	void (*print)(Model*, FILE*);
    void (*jsonize)(Model*, json_node*);
	double epsilon; // for finite differences
};

Model * new_Model( model_t type, const char *name, void *obj );

void free_Model( Model *model );

double Model_first_derivative( Model *model, Parameter* parameter, double eps );

void Model_first_derivatives( Model *model, Parameter* parameter, double eps, double* grad );

double Model_second_derivative( Model *model, Parameter* parameter, double* first, double eps );

double Model_mixed_derivative( Model *model, Parameter* p1, Parameter* p2 );

// Generic finite-difference Hessian over the flattened elements of `parameters`,
// in natural parameter space. Installed as the default Model->hessian. `out` must
// hold dim entries (HESSIAN_DIAGONAL) or dim*dim entries (HESSIAN_FULL).
void Model_hessian_fd( Model *model, const Parameters *parameters,
                       hessian_mode_t mode, double *out );

#pragma mark -
#pragma mark CatParameterModel

Model * new_CatParameterModel( const char* name, Parameters *parameters );

Parameter* new_CatParameter_from_json(json_node* node, Hashtable* hash);

#endif
