/*
 *  mg94.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/12/2015.
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

#ifndef mg94_h
#define mg94_h

#include <stdio.h>

#include "substmodel.h"

SubstitutionModel* new_MG94(Parameter* freqs, unsigned gen_code);

SubstitutionModel* new_MG94_with_values(Parameter* freqs, const double alpha,
                                        const double beta, const double kappa,
                                        unsigned gen_code);

SubstitutionModel* new_MG94_with_parameters(Parameter* freqs, Parameter* alpha,
                                            Parameter* beta, Parameter* kappa,
                                            unsigned gen_code);

#endif /* mg94_h */
