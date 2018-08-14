/*
 *  gtr.h
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

/**
 * @file gtr.h
 *
 * @brief The General Time-Reversible substitution model for nucleotides.
 *
 * This model is sometimes called the REV model, following Yang (1994), see references.
 * It was used in Lanave et al (1984), described in Tavare et al. 1886 and Rodriguez et al. 1990.
 * It is the most general reversible one, it has 6 substitution rates and 4 frequency
 * parameters.
 * We used the parametrization proposed by Yang (1994):
 * \f[
 * S = \begin{pmatrix}
 * \cdots & a & b & c \\
 * a & \cdots & d & e \\
 * b & d & \cdots & 1 \\
 * c & e & 1 & \cdots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\pi_A, \pi_C, \pi_G, \pi_T\right)
 * \f]
 * \f[
 * Q = \frac{S . \pi}{C} = \frac{1}{C} \begin{pmatrix}
 * -(a \pi_C + b \pi_G + c \pi_T) & a \pi_C & b \pi_G & c \pi_T \\
 * a \pi_A & -(a \pi_A + d \pi_G + e \pi_T) & d \pi_G & e \pi_T \\
 * b \pi_A & d \pi_C & -(b \pi_A + d \pi_C +\pi_T) & \pi_T \\
 * c \pi_A & e \pi_C & \pi_G & -(c \pi_A + e \pi_C + \pi_G) \\
 * \end{pmatrix}
 * \f]
 * \f[
 * C = \pi_A(a \pi_C + b \pi_G + c \pi_T) + \pi_C(a \pi_A + d \pi_G + e \pi_T) + \pi_G(b \pi_A + d \pi_C +\pi_T) + \pi_T(c \pi_A + e \pi_C + \pi_G)
 * \f]
 * Frequencies are reparametrized
 * \f{eqnarray*}
 * \phi_1 &=& \frac{\pi_A}{\pi_T}\\
 * \phi_2 &=& \frac{\pi_C}{\pi_T}\\
 * \phi_3 &=& \frac{\pi_G}{\pi_T}\\
 * \pi_T &=& \frac{1}{1+\sum_{i=1}^3{\phi_i}}
 * \f}
 * \f[
 * \pi = \left(\phi_1 \pi_T, \phi_2 \pi_T, \phi_3 \pi_T, \pi_T\right)
 * \f]
 * \f[
 * Q = \frac{S . \pi}{C}  = \frac{1}{C} \begin{pmatrix}
 * -(a \phi_2 \pi_T + b \phi_3 \pi_T + c \pi_T) & a \phi_2 \pi_T & b \phi_3 \pi_T & c \pi_T \\
 * a \phi_1 \pi_T & -(a \phi_1 \pi_T + d \phi_3 \pi_T + e \pi_T) & d \phi_3 \pi_T & e \pi_T \\
 * b \phi_1 \pi_T & d \phi_2 \pi_T & -(b \phi_1 \pi_T + d \phi_2 \pi_T +\pi_T) & \pi_T \\
 * c \phi_1 \pi_T & e \phi_2 \pi_T & \phi_3 \pi_T & -(c \phi_1 \pi_T + e \phi_2 \pi_T + \phi_3 \pi_T) \\
 * \end{pmatrix}
 * \f]
 * \f{eqnarray*}
 * C &=& \phi_1 \pi_T(a \phi_2 \pi_T + b \phi_3 \pi_T + c \pi_T) + \phi_2 \pi_T(a \phi_1 \pi_T + d \phi_3 \pi_T + e \pi_T) + \phi_3 \pi_T(b \phi_1 \pi_T + d \phi_2 \pi_T +\pi_T) + \pi_T(c \phi_1 \pi_T + e \phi_2 \pi_T + \phi_3 \pi_T)\\
 *   &=& \phi_1 \pi_T^2(a \phi_2 + b \phi_3 + c) + \phi_2 \pi_T^2(a \phi_1 + d \phi_3 + e) + \phi_3 \pi_T^2(b \phi_1 + d \phi_2  + 1) + \pi_T^2(c \phi_1 + e \phi_2 + \phi_3)\\
 *	 &=& \pi_T^2 (\phi_1(a \phi_2 + b \phi_3 + c) + \phi_2(a \phi_1 + d \phi_3 + e) + \phi_3 (b \phi_1 + d \phi_2  + 1) + (c \phi_1 + e \phi_2 + \phi_3))\\
 * \f}
 *
 *
 \f{eqnarray*}
 * \phi_1 &=& \frac{\pi_A}{\pi_T}\\
 * \phi_2 &=& \frac{\pi_C}{\pi_T}\\
 * \phi_3 &=& \frac{\pi_G}{\pi_T}\\
 * \pi_T &=& \frac{1}{1+\sum_{i=1}^3{\phi_i}}
 * \f}
 * \f[
 * \pi = \left(\phi_1 \pi_T, \phi_2 \pi_T, \phi_3 \pi_T, \pi_T\right)
 * \f]
 * \f[
 * Q = \frac{S . \pi}{C}  = \frac{1}{C} \begin{pmatrix}
 * -(a \phi_2 + b \phi_3 + c) & a \phi_2 & b \phi_3 & c \\
 * a \phi_1 & -(a \phi_1 + d \phi_3 + e) & d \phi_3 & e \\
 * b \phi_1 & d \phi_2 & -(b \phi_1 + d \phi_2 +1) & 1\\
 * c \phi_1 & e \phi_2 & \phi_3 & -(c \phi_1 + e \phi_2 + \phi_3) \\
 * \end{pmatrix}
 * \f]
 * \f{eqnarray*}
 * C &=& \pi_T (\phi_1(a \phi_2 + b \phi_3 + c) + \phi_2(a \phi_1 + d \phi_3 + e) + \phi_3 (b \phi_1 + d \phi_2  + 1) + (c \phi_1 + e \phi_2 + \phi_3))\\
 * \f}
 * \f{equation}
 * \frac{\partial \pi_T}{\partial \phi_1} = \frac{\partial \pi_T}{\partial \phi_2} = \frac{\partial \pi_T}{\partial \phi_3} = - \frac{1}{(1+\sum_{i=1}^3{\phi_i})^2} = -\pi_T^2
 * \f}
 * \f{eqnarray*}
 * \frac{\partial C}{\partial \phi_1} &=& 2 \pi_T (a \phi_2 + b \phi_3 + c)  - \pi_T C \\
 * \frac{\partial C}{\partial \phi_2} &=& 2 \pi_T (a \phi_1 + d \phi_3 + e)  - \pi_T C \\
 * \frac{\partial C}{\partial \phi_3} &=& 2 \pi_T (b \phi_1 + d \phi_2 + 1)  - \pi_T C
 * \f}
 *
 * \f{eqnarray*}
 * \frac{\partial}{\partial \phi_k}\frac{Q_{i,i}}{C} &=& -\frac{Q_{i,j}}{C^2}  \frac{\partial C}{\partial \phi_k}
 * \f}
 * \f{eqnarray*}
 * \frac{\partial}{\partial \phi_k}\frac{Q_{i,j}}{C} &=& \frac{C \frac{\partial Q_{i,j}}{\partial \phi_k} - Q_{i,j} \frac{\partial C}{\partial \phi_k}}{C^2}
 * \f}
 *
 */

#ifndef gtr_h
#define gtr_h

#include <stdio.h>

#include "substmodel.h"
/**
 * @brief A general time-reversible (GTR) substitution model.
 */

struct SubstitutionModel * new_GTR(Simplex* freqs);

struct SubstitutionModel * new_GTR_with_parameters( Simplex* freqs, Parameter* ac, Parameter* ag, Parameter* at, Parameter* cg, Parameter* ct, Parameter* gt);

struct SubstitutionModel * new_GTR_with_values( const double *freqs, const double *rates );

struct SubstitutionModel * new_GTR_with_simplexes( Simplex* freqs, Simplex* rates);

void gtr_update_Q( struct SubstitutionModel *m );

void gtr_simplexes_update_Q( struct SubstitutionModel *m );

void gtr_dQ(struct SubstitutionModel *m, int index, double* mat, double t);

#endif /* gtr_h */
