/*
 *  treelogio.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 29/8/12.
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
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "tree.h"
#include "filereader.h"
#include "mstring.h"
#include "treelogio.h"
#include "modelselection.h"

// TODO: Sould implement also the local clock models

int TreeLog_define_range( const char *filename, int *start, int *end, double alpha ){
	FileReader *reader = new_FileReader(filename,100);
	StringBuffer *buffer = new_StringBuffer(100);
	
	char *ptr = NULL;
	*start = 0;
	int nRate_current = 0;
	
	double IC_best_previous = INFINITY;
	double IC_best_current = INFINITY;
	
	double lnl_best_previous = -INFINITY;
	double lnl_best_current  = -INFINITY;
	
	int index = 0;
	int index_previous = 0;
	int index_current = 0;
	
	int nRate_start = 0;
	
	double p = 0; // for LRT
	
	bool discrete = false;
	
	while ( reader->read_line(reader) ) {
		ptr = reader->line;
		if ( String_start_with(reader->line, "tree", true) ) {
            
			int nRate = 0;
			if ( String_contains_str(ptr, "local=1") ) {
				nRate = String_contains_str_count(ptr, "local=1");
			}
			else if( String_contains_str(ptr, "class=") ){
				int pos = 0;
				while ( (pos = String_index_of_str(ptr, "class=") ) != -1 ) {
					ptr += pos + strlen("class=");
					StringBuffer_empty(buffer);
					while ( *ptr >= 48 && *ptr <= 57 ) {
						StringBuffer_append_char(buffer, *ptr);
						ptr++;
					}
					int c = atoi(buffer->c);
					if ( c > nRate) {
						nRate = c;
					}
					
				}
				discrete = true;
			}
			else {
				//error("BranchModel is neither local nor discrete\n");
			}
			
			nRate++;
			//fprintf(stderr, "%d nrate %d\n", index, nRate);
			
			
			if ( index == 0 ) {
				nRate_start = nRate;
			}
			
			
			StringBuffer_empty(buffer);
			ptr = reader->line;
			
			double IC = INFINITY;
			double lnl = INFINITY;
			
			//tree TREE1 [&LnL=-4154.174777,IC=8470.349553]]
            if ( String_contains_str(ptr, "AICc=") ) {
                ptr += String_index_of_str(ptr, "AICc=") + strlen("AICc=");
            }
            
			if ( String_contains_str(ptr, "IC=") ) {
				ptr += String_index_of_str(ptr, "IC=") + strlen("IC=");
			}
            while ( isspace(*ptr) ) ptr++;
            while ( (*ptr >= 48 && *ptr <= 57) || *ptr == '.' || *ptr == '-' ) {
                StringBuffer_append_char(buffer, *ptr);
                ptr++;
            }
            IC = atof(buffer->c);
			
			ptr = reader->line;
			StringBuffer_empty(buffer);
			if ( String_contains_str(ptr, "LnL=") ){
				ptr += String_index_of_str(ptr, "LnL=") + strlen("LnL=");
				while ( isspace(*ptr) ) ptr++;
				while ( (*ptr >= 48 && *ptr <= 57) || *ptr == '.' || *ptr == '-' ) {
					StringBuffer_append_char(buffer, *ptr);
					ptr++;
				}
				lnl = atof(buffer->c);
				
			}
            
			if( IC == INFINITY && lnl == INFINITY ){
				error("Could not read IC and LnL in TreeLog_define_range\n");
			}
            
			//fprintf(stdout, "%d %d - %f %f\n", nRate_previous, nRate_current, IC_best_previous, IC_best_current);
            
			// first batch of models
			if( nRate == nRate_start ){
				IC_best_previous = dmin(IC_best_previous, IC);
				lnl_best_previous = dmax(lnl_best_previous, lnl);
				index_previous = 0;
			}
			// second batch of models
			// we still have nothing to compare to
			else if ( nRate == nRate_start+1 ){
				if ( IC_best_current == INFINITY ) {
					index_current = index;
				}
				IC_best_current = dmin(IC_best_current, IC);
				lnl_best_current = dmax(lnl_best_current, lnl);
                
				nRate_current = nRate;
			}
			else if( nRate != nRate_current ) {
				//fprintf(stdout, "change model %d previous %d IC previous %f %f\n", index, index_previous, IC_best_current, IC_best_previous );
				if ( discrete ) {
					if ( IC_best_current < IC_best_previous ) {
						index_previous = index_current;
						IC_best_previous = IC_best_current;
						
						index_current = index;
						IC_best_current = IC;
						nRate_current = nRate;
					}
					else {
						fprintf(stderr, "Should not be reached unless we used the run was stopped (discrete)\n");
						break;
					}
				}
				else {
					p = LRT( lnl_best_previous, lnl_best_current, 0, 1);
					//fprintf(stderr, "LRT %f %f %f\n", lnl_best_previous, lnl_best_current, p);
					if ( p < alpha ) {
						index_previous = index_current;
						lnl_best_previous = lnl_best_current;
						
						index_current = index;
						lnl_best_current = lnl;
						nRate_current = nRate;
					}
					else {
						// no need to continue to read since the p < alpha
						break;
					}
				}
			}
			// same model
			else {
				IC_best_current  = dmin(IC_best_current, IC);
				lnl_best_current = dmax(lnl_best_current, lnl);
			}
			
			index++;
		}
	}
	free_FileReader(reader);
	free_StringBuffer(buffer);
	
	int numberOfRates = 0;
	// we have reached the hard limit
	// not tested
	if ( (discrete && IC_best_current < IC_best_previous) || (!discrete && p < alpha) ) {
		fprintf(stderr, "The very last set of models were better than the previous one (may have reached the hard limit during the GA search\n\n)");
		*start = index_current;
		*end = index;
		numberOfRates = nRate_current;
	}
	else {
		*start = index_previous;
		if ( index_current == 0 ) {
			*end = index;
		}
		else {
			*end = index_current-1;
		}
		numberOfRates = nRate_current-1;
	}
	return numberOfRates;
}

// Find the indexes of the tree that have the specified number of discrete clock
// assumes that trees with the number of rates are grouped together
void TreeLog_define_range_nclasses( const char *filename, int nclasses, int *start, int *end ){
	FileReader *reader = new_FileReader(filename,100);
	StringBuffer *buffer = new_StringBuffer(100);
	
	char *ptr = NULL;
    
    *start = -1;
    *end = -1;
	int index = 0;
    
	while ( reader->read_line(reader) ) {
		ptr = reader->line;
		if ( String_start_with(reader->line, "tree", true) ) {
            
			int nRate = 0;
			if ( String_contains_str(ptr, "local=1") ) {
				nRate = String_contains_str_count(ptr, "local=1");
			}
			else if( String_contains_str(ptr, "class=") ){
				int pos = 0;
				while ( (pos = String_index_of_str(ptr, "class=") ) != -1 ) {
					ptr += pos + strlen("class=");
					StringBuffer_empty(buffer);
					while ( *ptr >= 48 && *ptr <= 57 ) {
						StringBuffer_append_char(buffer, *ptr);
						ptr++;
					}
					int c = atoi(buffer->c);
					if ( c > nRate) {
						nRate = c;
					}
					
				}
			}
			else {
				error("BranchModel is neither local nor discrete\n");
			}
			
			nRate++;
			//fprintf(stderr, "%d nrate %d\n", index, nRate);
			
            if( nRate == nclasses && *start == -1 ){
                *start = index;
            }
            
            if( nRate != nclasses && *start != -1 ){
                *end = index-1;
                break;
            }
			
			index++;
		}
	}
	free_FileReader(reader);
	free_StringBuffer(buffer);
    
    if( *start != 1 && *end == -1 ){
        *end = index-1;
    }
}
