// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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

Sequences* new_Sequences_from_json(json_node* node, Hashtable* hash);

#endif
