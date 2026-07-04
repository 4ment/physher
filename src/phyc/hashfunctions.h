// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

/*
 **************************************************************************
 *                                                                        *
 *          General Purpose Hash Function Algorithms Library              *
 *                                                                        *
 * Author: Arash Partow - 2002                                            *
 * URL: http://www.partow.net                                             *
 * URL: http://www.partow.net/programming/hashfunctions/index.html        *
 *                                                                        *
 * Copyright notice:                                                      *
 * Free use of the General Purpose Hash Function Algorithms Library is    *
 * permitted under the guidelines and in accordance with the most current *
 * version of the Common Public License.                                  *
 * http://www.opensource.org/licenses/cpl1.0.php                          *
 *                                                                        *
 **************************************************************************
 */



#ifndef INCLUDE_HASHFUNCTIONS_H
#define INCLUDE_HASHFUNCTIONS_H


#include <stdio.h>


typedef unsigned int (*hash_function)(const void* data);

unsigned int JSHash2(const void* data);

unsigned int RSHash  ( const char* str, unsigned int len);
unsigned int JSHash  ( const char* str, unsigned int len);
unsigned int PJWHash ( const char* str, unsigned int len);
unsigned int ELFHash ( const char* str, unsigned int len);
unsigned int BKDRHash( const char* str, unsigned int len);
unsigned int SDBMHash( const char* str, unsigned int len);
unsigned int DJBHash ( const char* str, unsigned int len);
unsigned int DEKHash ( const char* str, unsigned int len);
unsigned int BPHash  ( const char* str, unsigned int len);
unsigned int FNVHash ( const char* str, unsigned int len);
unsigned int APHash  ( const char* str, unsigned int len);


#endif

