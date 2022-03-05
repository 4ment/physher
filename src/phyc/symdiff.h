//
//  symdiff.h
//  physher
//
//  Created by Mathieu Fourment on 2/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef symdiff_h
#define symdiff_h

#include <stdio.h>
#include <stdbool.h>

#include "parameters.h"

typedef struct ExpressionItem{
	char* expr;
	char* dev;
	char op;
	double value;
}ExpressionItem;

typedef struct ExpressionStack{
	ExpressionItem** expressions;
	size_t capacity;
	size_t count;
}ExpressionStack;

ExpressionStack* FillStack(const char* lpcsInput);

void free_Stack(ExpressionStack* stack);

char* DifferentiateStack(ExpressionStack* stack, int* nExpression, const char* dx);

char* differentiate(char* lpcsInput, const char* dx);

double differentiate2(ExpressionStack* stack, Parameters* list);

void Optimize(char* str);

#endif /* symdiff_h */
