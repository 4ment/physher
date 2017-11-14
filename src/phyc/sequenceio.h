/*
 *  sequenceio.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/2/10.
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

#ifndef _SEQUENCEIO_H_
#define _SEQUENCEIO_H_

#include "sequence.h"
#include "mstring.h"
#include "mjson.h"

Sequences * readSequences ( const char *infile );

Sequences * readFasta( const char *infile );

Sequences * readNexus ( const char *infile, int index );

Sequences * readPhylip ( const char *infile );

void Sequences_save_fasta( const Sequences *sequences, const char *filename );

void Sequences_save_nexus( const Sequences *sequences, const char *filename );

void Sequences_save_nexus_with_comment( const Sequences *sequences, const char *filename, char *comment );

void Sequences_save_phylip( const Sequences *sequences, const char *filename );


char * File_stringify( const char *filename );

StringBuffer * File_bufferize( StringBuffer *buffer, const char *filename );

Sequences* new_Sequences_from_json(json_node* node, Hashtable* hash);

#endif
