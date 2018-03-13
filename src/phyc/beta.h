//
//  beta.h
//  physher
//
//  Created by Mathieu Fourment on 15/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef beta_h
#define beta_h

#include <stdio.h>

double invbetai(double p, double a, double b);

double dbetaprime(double x, double alpha, double beta);

double dlogbetaprime(double x, double alpha, double beta);

#endif /* beta_h */
