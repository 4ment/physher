/*
 *  parameters.c
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

#include "parameters.h"

#include <assert.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "mstring.h"
#include "utils.h"
#include "parser.h"
#include "matrix.h"

#define STRINGIFY_PARAMETERS  "params"
#define STRINGIFY_VALUE "value"
#define STRINGIFY_CONSTRAINT "constraint"


struct _Constraint{
	double lower;
	double upper;
	bool lower_fixed;
	bool upper_fixed;
	double flower;
	double fupper;
};

struct _Parameters{
	Parameter **list;
	size_t count;
	size_t capacity;
};


#pragma mark -
#pragma mark Constraint

Constraint * new_Constraint( const double lower, const double upper ){
	Constraint *c = (Constraint *)malloc( sizeof(Constraint) );
	assert(c);
	c->lower = lower;
	c->upper = upper;
	c->lower_fixed = false;
	c->upper_fixed = false;
	
	c->flower = lower;
	c->fupper = upper;
	//( ( lower == DBL_MIN && upper == DBL_MAX ) ? true : false);
	return c;
}

void free_Constraint( Constraint *c ){
	free(c);
	c = NULL;
}

Constraint * clone_Constraint( Constraint *cnstr ){
	Constraint *new = new_Constraint(cnstr->lower, cnstr->upper);
	new->lower_fixed = cnstr->lower_fixed;
	new->upper_fixed = cnstr->upper_fixed;
	return new;
}

StringBuffer * Constraint_bufferize( StringBuffer *buffer, Constraint *c ){
	buffer = StringBuffer_append_format(buffer,"  (lower:%f)(upper:%f)", c->lower, c->upper);
	return buffer;
}

void * Constraint_SML_to_object( SMLNode node ){
	Constraint *c = NULL;
	SMLNode lower = NULL;
	SMLNode upper = NULL;
	if ( (lower = SML_get_element(node, "lower")) != NULL && (upper = SML_get_element(node, "upper")) != NULL ){
		c = new_Constraint( atof( SML_get_data(lower)), atof( SML_get_data(upper)) );
	}
	return c;
}

void Constraint_set_bounds( Constraint *constr, const double lower, const double upper ){
	if( constr != NULL ){
		constr->lower = lower;
		constr->upper = upper;
		constr->flower = lower;
		constr->fupper = upper;
	}
}

bool Constraint_lower_fixed( const Constraint *c ){
	return c->lower_fixed;
}

bool Constraint_upper_fixed( const Constraint *c ){
	return c->upper_fixed;
}

void Constraint_set_lower_fixed( Constraint *c, const bool fixed ){
	c->lower_fixed = fixed;
}


void Constraint_set_upper_fixed( Constraint *c, const bool fixed ){
	c->upper_fixed = fixed;
}


double Constraint_upper( const Constraint *c ){
	return c->upper;
}

double Constraint_lower( const Constraint *c ){
	return c->lower;
}

double Constraint_fupper( const Constraint *c ){
	return c->fupper;
}

double Constraint_flower( const Constraint *c ){
	return c->flower;
}

void Constraint_set_upper( Constraint *c, const double upper ){
	c->upper = upper;
}

void Constraint_set_lower( Constraint *c, const double lower ){
	c->lower = lower;
}

void Constraint_set_fupper( Constraint *c, const double fupper ){
	c->fupper = fupper;
}

void Constraint_set_flower( Constraint *c, const double flower ){
	c->flower = flower;
}

void compare_constraint( const Constraint *c1, const Constraint *c2 ){
	assert( c1->lower == c2->lower );
	assert( c1->upper == c2->upper );
}

#pragma mark -
#pragma mark Parameter

// Owns the constraint by default
Parameter * new_Parameter( const char *name, const double value, Constraint *constr ){
	return new_Parameter_with_postfix( name, "", value, constr );
}

Parameter * new_Parameter_with_postfix( const char *name, const char *postfix, const double value, Constraint *constr ){
	Parameter *p = (Parameter *)malloc( sizeof(Parameter) );
	assert(p);
	size_t name_len    = strlen(name);
	size_t postfix_len = strlen(postfix);
	
	p->name = (char*)malloc( sizeof(char)*(name_len+postfix_len+2) );
	assert(p->name);
	strcpy(p->name, name);
	if ( postfix_len > 0 ) {
		p->name[name_len] = '.';
		strcpy(p->name+name_len+1, postfix);
	}
	p->value = value;
	p->cnstr = constr;
	p->estimate = true;
	p->id = 0;
	p->listeners = new_ListenerList(1);
	p->refCount = 1;
	return p;
}

Parameter * clone_Parameter( Parameter *p ){
	Parameter *pnew = NULL;
	Constraint *cnstr = NULL;

	if( p->cnstr != NULL ){
		cnstr = clone_Constraint(p->cnstr);
	}
	
	pnew = new_Parameter( p->name, p->value,  cnstr );//pnew ref_count is reset to 1
	pnew->estimate = p->estimate;
	pnew->id = p->id;
	return pnew;
}

void free_Parameter( Parameter *p ){
	if(p->refCount == 1){
		free(p->name);
		p->listeners->free(p->listeners);
		if( p->cnstr != NULL ) free_Constraint(p->cnstr);
		free(p);
	}
	else{
		p->refCount--;
	}
}

char * Parameter_stringify( Parameter *p ){
	StringBuffer *buffer = new_StringBuffer( 100);
	
	buffer = Parameter_bufferize( buffer, p );
	
	free_StringBuffer( buffer );
	char *str = StringBuffer_tochar( buffer );
	return str;
}


StringBuffer * Parameter_bufferize( StringBuffer *buffer, Parameter *p ){
	char *pch = strrchr(p->name,'.');
	buffer = StringBuffer_append_string(buffer," (");
	
	if ( pch != NULL ){
		int l = pch - p->name;
		buffer = StringBuffer_append_range( buffer, p->name, 0,l);
	}
	else {
		buffer = StringBuffer_append_format(buffer,"%s", p->name);
	}
	
	buffer = StringBuffer_append_format(buffer,":\n  (%s:%f)\n", STRINGIFY_VALUE, p->value);	
	
	if( p->cnstr != NULL ){
		buffer = Constraint_bufferize( buffer, p->cnstr);
		buffer = StringBuffer_append_char(buffer,'\n');
	}
	
	buffer = StringBuffer_append_string(buffer," )");
	return buffer;
}

// write without the constraint, called when it is part of a set parameters with the same constraint
// chop the postfix if there is one
StringBuffer * Parameter_bufferize_parameter_set( StringBuffer *buffer, Parameter *p ){
	char *pch = strrchr(p->name,'.');
	if ( pch != NULL ){
		int l = pch - p->name;
		buffer = StringBuffer_append_string(buffer," (");
		buffer = StringBuffer_append_substring( buffer, p->name, l);
		buffer = StringBuffer_append_format(buffer,":%f)", p->value);
	}
	else {
		buffer = StringBuffer_append_format(buffer," (%s:%f)", p->name, p->value);
	}
	return buffer;
}

void * Parameter_SML_to_object( SMLNode node ){
	return Parameter_SML_to_object_with_postfix( node, "");
}

void * Parameter_SML_to_object_with_postfix( SMLNode node, const char *postfix ){
	Parameter *p = NULL;
	// shared constraint
	if ( SML_get_child_count( node ) == 0) {
		p = new_Parameter_with_postfix( SML_get_tag(node), postfix, atof(SML_get_data(node)), NULL);
	}
	else {
		Constraint *c = NULL;
		SMLNode lower = NULL;
		SMLNode upper = NULL;
		if ( (lower = SML_get_element(node, "lower")) != NULL ){
			upper = SML_get_element( node, "upper");
			c = new_Constraint( atof( SML_get_data(lower)), atof( SML_get_data(upper)) );
		}
		SMLNode value = SML_get_element( node, "value");
		p = new_Parameter_with_postfix( SML_get_tag(node), postfix, atof(SML_get_data(value)), c);
	}
	return p;
}

void Parameter_set_name( Parameter *p, const char *name ){
	if ( p->name != NULL ) {
		free(p->name);
	}
	p->name = String_clone(name);
}

char * Parameter_name( const Parameter *p ){
	return p->name;
}

void Parameter_set_value( Parameter *p, const double value ){
	p->value = value;
	p->listeners->fire(p->listeners, NULL, p->id);
}



double Parameter_value( const Parameter *p ){
	return p->value;
}

void assign_value( Parameter *p, const double value ){
	if ( p->cnstr == NULL ) 
		p->value = value;
	else{
		if ( value >= p->cnstr->lower &&  value <= p->cnstr->upper ) p->value = value;
		else if( value < p->cnstr->lower ) p->value = p->cnstr->lower;
		else if( value > p->cnstr->upper ) p->value = p->cnstr->upper;
		else {
			fprintf(stderr, "\n\n");
			Parameter_print(p);
			error("parameters.c:assign_value:should not be here\n");
		}
	}
#ifdef ASSIGN_DEBUG
	Parameter_die(p, value);
#endif
}

Constraint * Parameter_constraint( const Parameter *p ){
	return p->cnstr;
}

bool Parameter_estimate( const Parameter *p ){
	return p->estimate;
}

void Parameter_set_estimate( Parameter *p, const bool estimate ){
	p->estimate = estimate;
}

double Parameter_upper( const Parameter *p ){
	return Constraint_upper( p->cnstr );
}

double Parameter_lower( const Parameter *p ){
	return Constraint_lower( p->cnstr );	
}

void Parameter_set_upper( Parameter *p, const double value ){
	Constraint_set_upper( p->cnstr, value );
}

void Parameter_set_lower( Parameter *p, const double value ){
	Constraint_set_lower( p->cnstr, value );	
}

void Parameter_set_bounds( Parameter *p, const double lower, const double upper ){
	if ( p->cnstr == NULL ) {
		p->cnstr = new_Constraint(lower, upper);
	}
	else {
		Constraint_set_bounds( p->cnstr, lower, upper );
	}
}


void Parameter_warn( const Parameter *p, const double value ){
	if ( p->cnstr != NULL && (value < p->cnstr->lower || value > p->cnstr->upper ) ){
		fprintf(stderr, "Parameter out of bound\n%s\n",__FILE__);
		Parameter_print(p);
		fprintf(stderr, "\n");
	}
}

void Parameter_die( const Parameter *p, const double value ){
	if ( p->cnstr != NULL && (value < p->cnstr->lower || value > p->cnstr->upper ) ){
		fprintf(stderr, "Parameter out of bound\n%s\n",__FILE__);
		Parameter_print(p);
		fprintf(stderr, "\n");
		exit(1);
	}
}

double check_value( Constraint *cnstr, double value ){
	if ( cnstr == NULL ) {
		return value;
	}
	else{
		if ( value >= cnstr->lower &&  value <= cnstr->upper ) return value;
#ifdef ASSIGN_DEBUG
		fprintf(stderr, "Parameter out of bound\n%s\n",__FILE__);
		fprintf(stderr, "%f [%f - %f]\n", value, cnstr->lower, cnstr->upper);
		exit(0);
#endif
		if( value < cnstr->lower ) return cnstr->lower;
		else if( value > cnstr->upper ) return cnstr->upper;
		fprintf(stderr, "lower %f value = %f upper = %f\n", cnstr->lower ,value, cnstr->upper);
		error("parameters.c:check_value:should not be here\n");
	}
	return -1;
}



void Parameter_print( const Parameter *p ){
	fprintf(stderr, "%s %e [%e - %e] estimate = %d lower_fixed = %d upper_fixed = %d\n",p->name, p->value, p->cnstr->lower, p->cnstr->upper, p->estimate, p->cnstr->lower_fixed, p->cnstr->upper_fixed);
}

void compare_parameter( const Parameter *p1, const Parameter *p2 ){
	assert(p1->value == p2->value);
	assert( strcmp(p1->name, p2->name) == 0 );
	if( p1->cnstr != NULL && p2->cnstr != NULL )
		compare_constraint( p1->cnstr, p2->cnstr );
	else assert( p1->cnstr == NULL && p2->cnstr == NULL );
}

#pragma mark -
#pragma mark Parameters

Parameters * new_Parameters( const size_t capacity ){
	Parameters *ps = (Parameters *)malloc( sizeof(Parameters) );
	assert(ps);
	ps->list = (Parameter **)malloc( capacity * sizeof(Parameter *) );
	assert(ps->list);
	for (int i = 0; i < capacity; i++) {
		ps->list[i] = NULL;
	}
	ps->count = 0;
	ps->capacity = capacity;
	return ps;
}

void free_Parameters( Parameters *ps ){
	if(ps == NULL) return;
	
	for (int i = 0; i < ps->count; i++) {
		free_Parameter(ps->list[i]);
	}
	free(ps->list);
	free(ps);
	ps = NULL;
}

// The cloned parameters can have different capacities but same count of course
Parameters * clone_Parameters( Parameters *p ){
	Parameters * clone = new_Parameters( p->capacity );
	for (int i = 0; i < p->count; i++) {
		Parameters_move(clone, clone_Parameter( p->list[i] ) );
	}
	return clone;
}

StringBuffer * Parameters_bufferize( StringBuffer *buffer, Parameters *ps ){
    int i = 0;
    for ( i = 0; i < ps->count; i++) {
        StringBuffer_append_format(buffer, "%f%s", Parameters_value(ps, i),(ps->count-1 == i?"":","));
    }
    return buffer;
}

char * Parameters_stringify( Parameters *ps ){
	StringBuffer *buffer = new_StringBuffer( 100);
	
	buffer = Parameters_SML_bufferize( buffer, ps );
	
	free_StringBuffer( buffer );
	char *str = StringBuffer_tochar( buffer );
	return str;
}

StringBuffer * Parameters_SML_bufferize( StringBuffer *buffer, Parameters *ps ){
	int i = 0;
	
	char *postfix = ps->list[0]->name;
	char *pch = NULL;
	char *tempfix = NULL;
	for ( i = 0; i < ps->count; i++) {
		
		pch = strrchr( ps->list[i]->name,'.' );
		if ( pch != NULL ){
			tempfix = &ps->list[i]->name[pch - ps->list[i]->name+1];
			
			if ( strcmp( tempfix, postfix ) == 0 ){
				
			}
			else {
				if ( i!= 0) buffer = StringBuffer_append_string( buffer, ")\n"); 
				postfix = tempfix;
				buffer = StringBuffer_append_format( buffer, "(%s:\n", postfix);
				
			}
			buffer = Parameter_bufferize(buffer, ps->list[i]);
			buffer = StringBuffer_append_char(buffer, '\n');

		}
		else {
			buffer = Parameter_bufferize(buffer, ps->list[i]);
			if (ps->count-1 != i) buffer = StringBuffer_append_char(buffer, '\n');
		}
		
	}
	buffer = StringBuffer_append_string(buffer, ")\n");
	
	return buffer;
}


Parameters * Parameters_SML_to_object( SMLNode node ){
	SMLNode params = NULL;
	Parameters *ps =  NULL;
	int i = 0;
	
	char *name = SML_get_tag(node);
	
	// Set of parameters with the same constraint
	if ( (params = SML_get_element( node, "params")) != NULL ) {
		ps = new_Parameters(SML_get_child_count( params ));
		
		for ( i = 0; i < SML_get_child_count( params ); i++) {
			Parameters_add( ps, Parameter_SML_to_object( SML_get_element_by_index(params, i)) );
		}
		
	}
	else{
		ps = new_Parameters( SML_get_child_count( node ) );
		for ( i = 0; i < SML_get_child_count( node ); i++) {
			Parameters_add( ps, Parameter_SML_to_object_with_postfix( SML_get_element_by_index(node, i), name) );
		}
	}
	
	return ps;
}

Parameter * Parameters_at( const Parameters *p, const size_t index ){
	return p->list[index];
}

void Parameters_add( Parameters *ps, Parameter *p){
	if( ps->count == ps->capacity ){
		ps->list = realloc(ps->list, (ps->capacity+1) * sizeof(Parameter*) );
		assert(ps->list);
		ps->capacity++;
	}
	ps->list[ps->count++] = p;
	p->refCount++;
}

// add and do not increment refCount (move ownership)
void Parameters_move( Parameters *ps, Parameter *p){
	if( ps->count == ps->capacity ){
		ps->list = realloc(ps->list, (ps->capacity+1) * sizeof(Parameter*) );
		assert(ps->list);
		ps->capacity++;
	}
	ps->list[ps->count++] = p;
}

void Parameters_add_parameters(Parameters *dst, const Parameters *src){
	if(Parameters_count(src) == 0) return;
	if( dst->count+src->count > dst->capacity ){
		dst->list = realloc(dst->list, (dst->count+src->count) * sizeof(Parameter*) );
		assert(dst->list);
		dst->capacity = dst->count+src->count;
	}
	for (int i = 0; i < src->count; i++) {
		Parameters_add( dst, src->list[i] );
	}
}

void Parameters_add_free_parameters(Parameters *dst, const Parameters *src){
	if(Parameters_count(src) == 0) return;
	if( dst->count+src->count > dst->capacity ){
		dst->list = realloc(dst->list, (dst->count+src->count) * sizeof(Parameter*) );
		assert(dst->list);
		dst->capacity = dst->count+src->count;
	}
	for (int i = 0; i < src->count; i++) {
		if(src->list[i]->estimate){
			Parameters_add( dst, src->list[i] );
		}
	}
}

Parameters * pack_parameters(Parameters *ps, Parameter *p){
	ps->list = realloc(ps->list, ps->count * sizeof(Parameter*) );
	assert(ps->list);
	ps->capacity = ps->count;
	return ps;
}



void Parameters_set_name( Parameters *p, const size_t index, const char *name ){
	assert( index < p->count);
	Parameter_set_name(p->list[index], name);
}

char * Parameters_name( const Parameters *p, const size_t index ){
	return Parameter_name(p->list[index]);
}

size_t Parameters_count( const Parameters *p ){
    if( p == NULL ) return 0;
	return p->count;
}


size_t Parameters_capacity( const Parameters *p ){
	return p->capacity;
}

void Parameters_set_value( Parameters *p, const int index, const double value ){
	assert( index < p->count);
	Parameter_set_value(p->list[index], value);
}

void Parameters_set_all_value( Parameters *p, const double value ){
	for (int i = 0; i < Parameters_count(p); ++i) {
		Parameter_set_value(p->list[i], value);
	}
}

double Parameters_value( const Parameters *p, const int index ){
	return Parameter_value(p->list[index]);
}

bool Parameters_estimate( const Parameters *p, const int index ){
	return Parameter_estimate(p->list[index] );
}

void Parameters_set_estimate( Parameters *p, const bool estimate, const int index ){
	Parameter_set_estimate(p->list[index], estimate );
}

Constraint * Parameters_constraint( const Parameters *p, const int index ){
	return Parameter_constraint( Parameters_at(p, index) );
}

double Parameters_upper( const Parameters *p, const int index ){
	return Parameter_upper( Parameters_at(p, index) );
}

double Parameters_lower( const Parameters *p, const int index ){
	return Parameter_lower( Parameters_at(p, index) );	
}

void Parameters_set_upper( Parameters *p, const int index, const double value ){
	Parameter_set_upper( Parameters_at(p, index), value );
}

void Parameters_set_lower( Parameters *p, const int index, const double value ){
	Parameter_set_lower( Parameters_at(p, index), value );	
}

void Parameters_set_bounds( Parameters *p, const int index, const double lower, const double upper ){
	Parameter_set_bounds( Parameters_at(p, index), lower, upper );
}


bool Parameters_is_at_boundry( const Parameters *p, const int index, double precision ){
    return ( Parameter_value(p->list[index]) - Parameters_lower(p, index) < precision
            || Parameters_upper(p, index) - Parameter_value(p->list[index]) < precision );
    
}

// not cloned
Parameters * get_sub_parameters( Parameters *p, const int start, const int end ){
	Parameters * sub = new_Parameters( end - start + 1 );
	for (int i = start; i <= end; i++) {
		Parameters_add(sub, p->list[i]);
	}
	return sub;
}

void Parameters_print( Parameters *ps ){
	fprintf(stderr, "==================================\n");
	for (int i = 0; i < ps->count; i++) {
		Parameter_print(ps->list[i]);
	}
	fprintf(stderr, "----------------------------------\n");
}

// The parameters are not cloned
// assume that .height are consecutive, followed by .bl...
Parameters * Parameters_optimizable( Parameters *p, int **map, const char postfix[] ){
	if( p == NULL ) return NULL;
	
	int count = Parameters_count_optimizable(p, postfix);
	
	if( count == 0 ) return NULL;

	*map = malloc( sizeof(int) * count );
	assert(map);
	
	Parameters * sub = new_Parameters( count );
	count = 0;

	int index = 0; 

	for ( int i = 0; i < p->count; i++) {
		if ( postfix != NULL ) {
			if( index == 0 && strstr(p->list[i]->name, postfix) != NULL )
				index = i;
		}
			
		
		if( p->list[i]->estimate ){
			if ( postfix != NULL ){
				if(strstr(p->list[i]->name, postfix) == NULL ) continue;
			}
			
			Parameters_add(sub, p->list[i]);
			(*map)[count++] = i - index;
		}
	}
	return sub;
}

void Parameters_store_value( const Parameters *p, double *store ){
	for (int i = 0; i < Parameters_count(p); i++) {
		store[i] = Parameters_value(p, i);
	}
}

void Parameters_restore_value( Parameters *p, const double *store ){
	for (int i = 0; i < Parameters_count(p); i++) {
		Parameters_set_value(p, i, store[i]);
	}
}


int Parameters_count_optimizable( const Parameters *p, const char postfix[] ){
	if( p == NULL ) return 0;
	int count = 0;
	for ( int i = 0; i < p->count; i++) {
		if( p->list[i]->estimate ){
			if ( postfix != NULL ){
				if(strstr(p->list[i]->name, postfix) != NULL ) count++;
			}
			else count++;
		}
	}
	return count;
}

void Parameters_remove( Parameters *params, size_t index ){
	assert(index < params->count);
	if( params->count !=0 ){
		free_Parameter( params->list[index] );
        
        for( size_t i = index+1; i < params->count; i++ ){
            params->list[i-1] = params->list[i];
        }
		params->list[params->count-1] = NULL;
		params->count--;
	}
}

void Parameters_pop( Parameters *params ){
	if( params->count !=0 ){
		free_Parameter( params->list[params->count-1] );
		params->list[params->count-1] = NULL;
		params->count--;
	}
}

void compare_parameters( const Parameters *p1, const Parameters *p2 ){
	assert( p1->count    == p2->count );
	assert( p1->capacity == p2->capacity );
	for (int i = 0; i < p1->count; i++) {
		compare_parameter( p1->list[i], p2->list[i]);
	}
}

void Parameters_swap( Parameter **a, Parameter **b ){
	Parameter *p = *a;
	*a = *b;
	*b = p;
}

void Parameters_swap_index( Parameters *ps, unsigned a, unsigned b ){
	Parameter *p = ps->list[a];
	ps->list[a] = ps->list[b];
	ps->list[b] = p;
}

void Parameters_sort_from_ivector( Parameters *p, int *s ){
	bool done = false;
	size_t size = p->count;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( s[i] > s[i+1] ) {
				done = false;
				swap_int( &s[i], &s[i+1] );
				Parameters_swap_index(p, i, i+1);
			}
		}
		size--;
	}
}

#pragma mark -
#pragma mark DiscreteParameter

static void _free_DiscreteParameter( DiscreteParameter *p ){
	if(p->refCount == 1){
		free(p->name);
		p->listeners->free(p->listeners);
		free(p);
	}
	else{
		p->refCount--;
	}
}

static DiscreteParameter * _clone_DiscreteParameter( DiscreteParameter *p ){
	DiscreteParameter *pnew = new_DiscreteParameter( p->name, p->length);
	pnew->id = p->id;
	memcpy(pnew->values, p->values, sizeof(unsigned)*p->length);
	return pnew;
}

static void _set_value_discrete(DiscreteParameter* p, int index, unsigned value){
	assert(index < p->length);
	p->values[index] = value;
	p->listeners->fire(p->listeners, NULL, index);
}

static void _set_values_discrete(DiscreteParameter* p, const unsigned* values){
	memcpy(p->values, values, p->length);
	p->listeners->fire(p->listeners, NULL, -1);
}

DiscreteParameter * new_DiscreteParameter( const char *name, int dim ){
	return new_DiscreteParameter_with_postfix(name, "", dim);
}

DiscreteParameter * new_DiscreteParameter_with_postfix( const char *name, const char *postfix, int dim ){
	return new_DiscreteParameter_with_postfix_values(name, postfix, NULL, dim);
}

DiscreteParameter * new_DiscreteParameter_with_values( const char *name, const unsigned* values, int dim ){
	return new_DiscreteParameter_with_postfix_values(name, "", values, dim);
}

DiscreteParameter * new_DiscreteParameter_with_postfix_values( const char *name, const char *postfix, const unsigned* values, int dim ){
	DiscreteParameter *p = (DiscreteParameter *)malloc( sizeof(DiscreteParameter) );
	assert(p);
	size_t name_len    = strlen(name);
	size_t postfix_len = strlen(postfix);
	
	p->name = (char*)malloc( sizeof(char)*(name_len+postfix_len+2) );
	assert(p->name);
	strcpy(p->name, name);
	if ( postfix_len > 0 ) {
		p->name[name_len] = '.';
		strcpy(p->name+name_len+1, postfix);
	}
	if(values != NULL){
		p->values = clone_uivector(values, dim);
	}
	else{
		p->values = uivector(dim);
	}
	p->set_value = _set_value_discrete;
	p->set_values = _set_values_discrete;
	p->free = _free_DiscreteParameter;
	p->clone = _clone_DiscreteParameter;
	p->id = 0;
	p->listeners = new_ListenerList(1);
	p->refCount = 1;
	return p;
}

static void dummyUpdate( Model *self, Model *model, int index ){}
static double _logP(Model *model){return 0;}
static double _dlogP(Model *model, const Parameter* p){return 0;}

#pragma mark -


Model * new_Model( const char *name, void *obj ){
	Model *model = (Model*)malloc(sizeof(Model));
	assert(model);
	model->name = String_clone(name);
	model->obj = obj;
	model->logP = _logP;
	model->dlogP = _dlogP;
	model->update  = dummyUpdate;
	model->free = free_Model;
	model->data = NULL;
	model->listeners = new_ListenerList(1);
	model->clone = NULL;
	model->get_free_parameters = NULL;
	model->ref_count = 1;
	return model;
}

void free_Model( Model *model ){
	assert(model->ref_count >= 1);
	if(model->ref_count == 1){
		free(model->name);
		model->listeners->free(model->listeners);
		free(model);
	}
	else{
		model->ref_count--;
	}
	
}

#pragma mark -

static void _free_ListenerList( ListenerList *listeners ){
	free(listeners->models);
	free(listeners);
}

static void _ListenerList_fire( ListenerList *listeners, Model* model, int index){
	for ( int i = 0; i < listeners->count; i++ ) {
		listeners->models[i]->update( listeners->models[i], model, index );
	}
}

static void _ListenerList_remove( ListenerList *listeners, Model* model ){
	int i = 0;
	for ( ; i < listeners->count; i++ ) {
		if ( listeners->models[i] == model ) {
			break;
		}
	}
	if ( i == listeners->count) {
		return;
	}
	i++;
	for ( ; i < listeners->count; i++ ) {
		listeners->models[i-1] = listeners->models[i];
	}
	listeners->models[listeners->count-1] = NULL;
	listeners->count--;
}

static void _ListenerList_remove_all( ListenerList *listeners ){
	int i = 0;
	for ( ; i < listeners->count; i++ ) {
		listeners->models[i] = NULL;
	}
	listeners->count = 0;
}

static void _ListenerList_add( ListenerList *listeners, Model *model ){
	if ( listeners->count == listeners->capacity) {
		listeners->capacity++;
		listeners->models = realloc(listeners->models, listeners->capacity*sizeof(Model*));
	}
	listeners->models[listeners->count] = model;
	listeners->count++;
}
ListenerList * new_ListenerList( const unsigned capacity ){
	ListenerList *listeners = (ListenerList*)malloc( sizeof(ListenerList));
	assert(listeners);
	listeners->capacity = (capacity == 0 ? 1 : capacity);
	listeners->count = 0;
	listeners->models = (Model**)malloc( listeners->capacity * sizeof(Model*));
	assert(listeners->models);
	listeners->free = _free_ListenerList;
	listeners->add = _ListenerList_add;
	listeners->remove = _ListenerList_remove;
	listeners->removeAll = _ListenerList_remove_all;
	listeners->fire = _ListenerList_fire;
	return listeners;
}



