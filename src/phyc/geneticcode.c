// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "geneticcode.h"

#include <stdio.h>

// if we want to print codon should end woith \0
void GenticCode_encoding_to_codon_string( int encoding, int genetic_code, char *codon ){
    if( encoding >= NUMBER_OF_CODONS[genetic_code] ){
        codon[0] = '-';
        codon[1] = '-';
        codon[2] = '-';
    }
    int count = 0;
    int i = 0;
    for ( ; i < 64; i++ ) {
        if( GENETIC_CODE_TABLES[genetic_code][i] == '*' ) continue;
        if(count == encoding ){
            break;
        }
        count++;
    }
    
    codon[0] = CODON_TRIPLETS[i][0];
    codon[1] = CODON_TRIPLETS[i][1];
    codon[2] = CODON_TRIPLETS[i][2];
}
