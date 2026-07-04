// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef utilsio_h
#define utilsio_h

#include <stdio.h>

#include "matrix.h"

Vector* read_log_column_with_id( const char *filename, size_t burnin, const char* id );

Vector** read_log_column_with_ids( const char *filename, size_t burnin, const char** tags, size_t tag_count );

#endif /* utilsio_h */
