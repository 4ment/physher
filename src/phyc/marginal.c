//
//  marginal.c
//  physher
//
//  Created by Mathieu Fourment on 18/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "marginal.h"

#include <float.h>

#include "utils.h"
#include "matrix.h"
#include "descriptivestats.h"
#include "statistics.h"

double log_harmonic_mean(const Vector** vecvalues){
	double sum = 0;
	double* values = Vector_data(vecvalues);
	size_t size = Vector_length(vecvalues);
	for(int i = 0; i < size; i++){
		sum += values[i];
	}
	
	double denominator = -DBL_MAX;
		
	for(int i = 0; i < size; i++){
		denominator = logaddexp(denominator, sum - values[i]);
	}
	return sum - denominator + log(size);
}

// temperatures should be sorted in increasing order
double log_marginal_stepping_stone(const Vector** values, size_t temp_count, const double* temperatures, double* lrssk){
	double lrss = 0;
	for(int i = 1; i < temp_count; i++){
		double tempdiff = temperatures[i]-temperatures[i-1];
		double logmaxll = dmax_vector(Vector_data(values[i]), Vector_length(values[i]));
		double* previousll = Vector_data(values[i-1]);
		double temp = -DBL_MAX;
		int size = Vector_length(values[i-1]);
		for (int j = 0; j < size; j++) {
			temp = logaddexp(temp, tempdiff*(previousll[j]-logmaxll));
		}
		lrssk[i] = tempdiff*logmaxll + temp - log(size);
		lrss += lrssk[i];
	}
	return lrss;
}

// temperatures should be sorted in increasing order
double log_marginal_path_sampling(const Vector** values, size_t temp_count, const double* temperatures, double* lrpsk){
	double lrps = 0;
	double* means = dvector(temp_count);
	for(int i = 0; i < temp_count; i++){
		means[i] = dmean(Vector_data(values[i]), Vector_length(values[i]));
	}
	
	for(int i = 1; i < temp_count; i++){
		double weight = temperatures[i]-temperatures[i-1];
		lrpsk[i] = weight * (means[i]+means[i-1])/2.0;
		lrps += lrpsk[i];
	}
	free(means);
	return lrps;
}

double log_marginal_path_sampling_modified(const Vector** values, size_t temp_count, const double* temperatures, double* lrpsk){
	double lrps = 0;
	double* means = dvector(temp_count);
	double* variances = dvector(temp_count);
	for(int i = 0; i < temp_count; i++){
		means[i] = dmean(Vector_data(values[i]), Vector_length(values[i]));
		variances[i] = variance(Vector_data(values[i]), Vector_length(values[i]), means[i]);
	}
	
	for(int i = 1; i < temp_count; i++){
		double weight = temperatures[i]-temperatures[i-1];
		lrpsk[i] = weight * (means[i]+means[i-1])/2.0 - ((weight*weight)/12)*(variances[i]-variances[i-1]);
		lrps += lrpsk[i];
	}
	free(means);
	free(variances);
	return lrps;
}
