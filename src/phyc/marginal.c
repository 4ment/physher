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
		double temp = 0;
		int size = Vector_length(values[i-1]);
		for (int j = 0; j < size; j++) {
			temp += exp(tempdiff*(previousll[j]-logmaxll));
		}
		lrssk[i] = tempdiff*logmaxll + log(temp/size);
		lrss += lrssk[i];
	}
	return lrss;
}
