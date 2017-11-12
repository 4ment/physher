/*
 *  filereader.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 12/14/11.
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

#endif
