/*
 *  sa.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 4/7/12.
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

#include "sa.h"

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "matrix.h"
#include "random.h"

#define BOLTZMANN_CONSTANT (1.3806503e-23)

SA * new_SA( const unsigned int state_size ){
	SA *sa = (SA*)malloc( sizeof(SA) );
	assert(sa);
	
	sa->temperature = 0.0;
	sa->initial_temp = 1e-10;
	sa->final_temp = 0.0;
	sa->temp_step = 0.0;
	sa->temp_freq = 0;
	
	sa->state_size = state_size;
	
	sa->state = new_State(state_size);
	sa->proposed_state = new_State(state_size);
	
	sa->max_iterations = 100;
	sa->iteration = 0;
	
	sa->tol = 0.001;
	sa->use_tol = false;
	
	sa->max_energy = 0;
	sa->best_state = uivector(state_size);
	
	sa->use_max_no_improvement = false;
	sa->max_no_improvement = 20;
	sa->n_without_improvement = 0;
	
	sa->mutate = NULL;
	sa->evaluate = NULL;
	sa->accept = NULL;
	
	sa->termination = sa_chatty_termination;
	
	
	sa->log = NULL;
	
	sa->output   = NULL;
	sa->filename = NULL;
	sa->dump = false;
	sa->dump_frequency = 0;
	
	sa->lookup = new_Hashtable_string(100);
	
	sa->data = NULL;
	
	sa->initialized = false;
	
	sa->start_time = 0;
	sa->current_time = 0;
	
	sa->eval = 0;

	
	return sa;
}


void free_SA( SA *sa ){
	free_State(sa->state);
	free_State(sa->proposed_state);

	free_Hashtable( sa->lookup);
	free(sa->filename);
	free(sa->best_state);
	free(sa);
}

State * new_State( const unsigned state_size ){
	State *state = (State*)malloc(sizeof(State));
	assert(state);
	state->encoding = uivector(state_size);
	state->size = state_size;
	state->data = NULL;
	state->energy = 0;
	return state;
}

void free_State( State *state ){
	free(state->encoding);
	free(state);
}

State * clone_State( const State *state ){
	State *newState = new_State(state->size);
	memcpy(newState->encoding, state->encoding, state->size*sizeof(unsigned*));
	newState->data = state->data;
	newState->energy = state->energy;
	return newState;
}



void sa( SA *sa ){
	State *temp_state = NULL;
	
	while( !sa->termination(sa) ){
		sa->iteration++;
		
		if (sa->temp_freq == -1) {
			sa->temperature = sa->initial_temp 
				+ ((double)sa->iteration/sa->max_iterations)
				* (sa->final_temp - sa->initial_temp);
		}
		else {
			if ( sa->temperature > sa->final_temp && sa->iteration % sa->temp_freq == 0 ) {
				sa->temperature -= sa->temp_step;
			}
		}
		
		sa->mutate(sa, sa->state, sa->proposed_state);
		sa->evaluate(sa, sa->proposed_state);
		
		if( sa->accept(sa, sa->state, sa->proposed_state) ){
			temp_state = sa->state;
			sa->state = sa->proposed_state;
			sa->proposed_state = temp_state;
		}
	}
}

bool ga_sa_boltzmann_acceptance( SA	*sa, State *state, State *proposed_state ) {
	if( proposed_state->energy < state->energy ){
		return false;
	}
	double p = exp( (proposed_state->energy - state->energy )/(BOLTZMANN_CONSTANT*sa->temperature));
	return random_double() < p;
}

bool sa_quiet_termination( SA *sa ){
	if( sa->max_iterations <= sa->iteration) return true;
	if( sa->use_max_no_improvement && sa->max_no_improvement == sa->n_without_improvement ) return true;
	return false;
}


bool sa_chatty_termination( SA *sa ){
	if (sa->iteration%1 == 0){
		time(&sa->current_time);
		double t = difftime(sa->current_time, sa->start_time);
		fprintf(stderr, "-----------------------------------------------------------------\n");
		fprintf(stderr, "Iteration = %d\n", sa->iteration);
		fprintf(stderr, "Number of iteration without improvement = %d\n", sa->n_without_improvement);		
		fprintf(stderr, "Max energy = %f\n", sa->max_energy);
		fprintf(stderr, "Total runtime ");
		print_pretty_time(stderr,t);
		fprintf(stderr, "Number of unique models evaluated %u (total %u)\n\n", Hashtable_length(sa->lookup), sa->eval);
    }
	
	if( sa->max_iterations <= sa->iteration) return true;
	if( sa->use_max_no_improvement && sa->max_no_improvement == sa->n_without_improvement ) return true;
	
	return ( sa->temperature < sa->final_temp);
		
	
	return false;
}
