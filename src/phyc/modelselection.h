/*
 *  ic.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/22/11.
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

#ifndef _MODEL_SELECTION_H_
#define _MODEL_SELECTION_H_

typedef enum information_criterion{
	INFORMATION_CRITERION_AIC,
	INFORMATION_CRITERION_AICc,
	INFORMATION_CRITERION_BIC,
	INFORMATION_CRITERION_HQ
}information_criterion;

static const char * const INFORMATION_CRITERION[4] = {
	"AIC",
	"AICc",
	"BIC",
	"HQ"
};


double AIC( const double lk, const int k );

double AICc( const double lk, const int k, const int n );

double BIC( const double lk, const int k, const int n );

double HQ( const double lk, const int k, const int n );

double LRT( const double lk0, const double lk1, const unsigned n0, const unsigned n1  );

#endif
