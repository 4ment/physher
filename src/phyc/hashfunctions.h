/*
 *  hashfunctions.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/5/10.
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

