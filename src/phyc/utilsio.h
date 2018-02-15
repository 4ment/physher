//
//  utilsio.h
//  physher
//
//  Created by Mathieu Fourment on 16/02/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef utilsio_h
#define utilsio_h

#include <stdio.h>

#include "matrix.h"

Vector* read_log_column_with_id( const char *filename, size_t burnin, const char* id );

Vector** read_log_column_with_ids( const char *filename, size_t burnin, const char** tags, size_t tag_count );

#endif /* utilsio_h */
