//
//  parametersio.h
//  physher
//
//  Created by Mathieu Fourment on 6/03/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef parametersio_h
#define parametersio_h

#include <stdio.h>
#include "parameters.h"
#include "matrix.h"

Vector** read_log_for_parameters( const char *filename, size_t burnin, size_t* count, Parameters* params );

Vector** read_log_for_parameters_t( const char *filename, size_t burnin, Parameters* params );

Vector** read_log_for_names_t( const char *filename, size_t burnin, char** params, size_t paramCount );

#endif /* parametersio_h */
