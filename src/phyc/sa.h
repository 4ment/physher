// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef Math_sa_h
#define Math_sa_h


#include <time.h>

#include "sa.h"
#include "utils.h"
#include "hashtable.h"



typedef struct _SA SA;

typedef struct State{
	unsigned int *encoding;
	unsigned int size;
	double energy;
	
	void *data;
}State;

typedef bool (*sa_accept)( SA *sa, State *state, State *proposed_state );
typedef void (*sa_mutate)( SA *sa, State *state, State *proposed_state );
typedef double (*sa_evaluate)( SA *sa, State *state );
typedef void (*sa_log)( SA *sa );
typedef bool (*sa_termination)( SA *sa );

struct _SA{
	double temperature;
	double initial_temp;
	double final_temp;
	double temp_step;
	int temp_freq;
	
	double cooling;
	
	unsigned int state_size;
	
	State *state;
	State *proposed_state;
	
	unsigned int max_iterations;
	unsigned int iteration;
	
	double tol;
	bool use_tol;
	
	double max_energy;
	unsigned int *best_state;
	
	bool use_max_no_improvement;
	int max_no_improvement;
	int n_without_improvement;
	
	
	sa_accept accept;
	sa_mutate mutate;
	sa_evaluate evaluate;
	
	sa_termination termination;
	
	sa_log log;
	
	FILE *output;
	char *filename;
	bool dump;
	int dump_frequency;
	
	Hashtable *lookup;
	
	void *data;
	
	bool initialized;
	
	time_t start_time;
	time_t current_time;
	
	unsigned eval;
};

SA * new_SA( const unsigned int state_size );

void free_SA( SA *sa );

State * new_State( const unsigned state_size );

void free_State( State *state );

State * clone_State( const State *state );
	


bool sa_chatty_termination( SA *sa );
bool sa_quiet_termination( SA *sa );

#endif
