// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef parametersio_h
#define parametersio_h

#include <stdio.h>
#include "parameters.h"
#include "matrix.h"

Vector** read_log_for_parameters( const char *filename, size_t burnin, size_t* count, Parameters* params );

Vector** read_log_for_parameters_t( const char *filename, size_t burnin, Parameters* params );

Vector** read_log_for_names_t( const char *filename, size_t burnin, char** params, size_t paramCount );

Vector** read_log_for_parameter_t( const char *filename, size_t burnin, Parameter* parameter );

#endif /* parametersio_h */
