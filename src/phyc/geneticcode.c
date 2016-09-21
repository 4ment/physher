/*
 *  geneticcode.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 13/03/2014.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include <stdio.h>


#include "geneticcode.h"

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
