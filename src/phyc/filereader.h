// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _FILE_READER_H_
#define _FILE_READER_H_

#include <stdio.h>

#include "mstring.h"

struct _FileReader;
typedef struct _FileReader FileReader;

struct _FileReader {
	FILE *file;
	StringBuffer *buffer;
	char *line;
    char b[513];
    size_t index;
    
    bool (*read_line)( FileReader * );
};

FileReader * new_FileReader( const char *filename, const unsigned buffer_size );

FileReader * new_FileReader_with_mode( const char *filename, const unsigned buffer_size, const char *mode );

void free_FileReader( FileReader *reader );

double ** FileReader_csv_double( const char *filename, int nrow, int ncol );

char* load_file(const char *filename);

// Read a whole stream (e.g. stdin) into a NUL terminated string.
char* load_stream(FILE *file);

#endif
