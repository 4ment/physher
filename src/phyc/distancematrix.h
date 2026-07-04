// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef __PhyC__distancematrix__
#define __PhyC__distancematrix__

#include <stdio.h>

#include "sequence.h"
#include "sitepattern.h"
#include "matrix.h"

typedef enum distancematrix_model{
    DISTANCE_MATRIX_UNCORRECTED,
    DISTANCE_MATRIX_JC69,
    DISTANCE_MATRIX_K2P,
    
    DISTANCE_MATRIX_KIMURA
}distancematrix_model;

double ** Sequences_distance( const Sequences *sequences, distancematrix_model model );

float ** Sequences_distance_float( const Sequences *sequences, distancematrix_model model );

double ** SitePattern_distance( const SitePattern *patterns, distancematrix_model model );

Matrix* create_DistanceMatrix_from_json( json_node* node, Hashtable* hash );

#endif /* defined(__PhyC__distancematrix__) */
