/*
 *  gensubst.c
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

#include "gensubst.h"

#include "matrix.h"


bool isSymmetric(const unsigned *matrix, size_t dim){
    for ( int i = 0; i < dim; i++ ) {
        for ( int j = i+1; j < dim; j++ ) {
            if( matrix[i*dim+j] != matrix[j*dim+i]){
                return false;
            }
        }
    }
    return true;
}


void _nonreversible_update_Q( SubstitutionModel *m ){
	const double* freqs = m->simplex->get_values(m->simplex);
    int index = 0;
	const unsigned* model = m->model->values;
    for ( int i = 0; i < m->nstate; i++ )  {
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++]) * freqs[j];
        }
        j++;
        for ( ; j < m->nstate; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++]) * freqs[j];
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
    if(m->normalize)normalize_Q( m->Q, freqs, m->nstate );
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    
    m->need_update = false;
}

void _nonreversible_update_Q2( SubstitutionModel *m ){
	const double* freqs = m->simplex->get_values(m->simplex);
    int index = 0;
	const unsigned* model = m->model->values;
    for ( int i = 0; i < m->nstate; i++ )  {
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++]) * freqs[j];
        }
        j++;
        for ( ; j < m->nstate; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++]) * freqs[j];
        }
    }
    
    for ( int i = 0; i < m->nstate-1; i++ )  {
        m->Q[m->nstate-1][i] = 0;
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[m->nstate-1][i] += m->Q[i][j];
        }
        for ( j++ ; j < m->nstate; j++ ) {
            m->Q[m->nstate-1][i] += m->Q[i][j];
        }
        for ( j = 0; j < m->nstate-1; j++ ) {
            if( i != j ){
                m->Q[m->nstate-1][i] -= m->Q[j][i];
            }
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
    if(m->normalize)normalize_Q( m->Q, freqs, m->nstate );
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    
    m->need_update = false;
}


void _reversible_update_Q( SubstitutionModel *m ){
	const unsigned* model = m->model->values;
	double temp;
	int index = 0;
	
	if(m->simplex != NULL){
		const double* freqs = m->simplex->get_values(m->simplex);
		for ( int i = 0; i < m->nstate; i++ )  {
			index += i+1;
			for ( int j = i + 1; j < m->nstate; j++ ) {
				if(m->relativeTo != model[index]){
					temp = Parameters_value(m->rates, model[index]);
				}
				else{
					temp = 1;
				}
				index++;
				m->Q[i][j] = temp * freqs[j];
				m->Q[j][i] = temp * freqs[i];
			}
		}
		make_zero_rows( m->Q, m->nstate);
		if(m->normalize)
			normalize_Q( m->Q, freqs, m->nstate );
	}
	else{
		for ( int i = 0; i < m->nstate; i++ )  {
			index += i+1;
			for ( int j = i + 1; j < m->nstate; j++ ) {
				if(m->relativeTo != model[index]){
					temp = Parameters_value(m->rates, model[index]);
				}
				else{
					temp = 1;
				}
				index++;
				m->Q[i][j] = temp;
				m->Q[j][i] = temp;
			}
		}
		make_zero_rows( m->Q, m->nstate);
	}
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    
    m->need_update = false;
}

SubstitutionModel * new_GeneralModel_with_parameters( DiscreteParameter* model, const Parameters* rates, Simplex* freqs, int relativeTo, bool normalize ){
    //bool sym = isSymmetric(model->values, freqs->K);
		size_t dim = sqrt(model->length);
	bool sym = isSymmetric(model->values, dim);
    SubstitutionModel *m = NULL;
    if(sym){
        m = create_general_model("GENERAL", REVERSIBLE, freqs, dim);
        m->update_Q = _reversible_update_Q;
    }
	else{
        m = create_general_model("UREVGENERAL", NONREVERSIBLE, freqs, dim);
        m->update_Q = _nonreversible_update_Q;
        m->reversible = false;
    }
    m->normalize = normalize;
	m->relativeTo = relativeTo;
    
	m->model = model;
	model->refCount++;
	
	m->rates = new_Parameters(Parameters_count(rates));
	Parameters_add_parameters(m->rates, rates);
    
    return m;
}

#pragma mark -
#pragma mark Reversible


static void _rev_update_Q( SubstitutionModel *m );
static void _rev_relative_update_Q( SubstitutionModel *m );

void _rev_update_Q( SubstitutionModel *m ){
    int index = 0;
    double temp;
	const unsigned* model = m->model->values;
    for ( int i = 0; i < m->nstate; i++ )  {
        for ( int j = i + 1; j < m->nstate; j++ ) {
            temp = Parameters_value(m->rates, model[index++]);
            m->Q[i][j] = temp;// * freqs[j];
            m->Q[j][i] = temp;// * freqs[i];
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
    //normalize_Q( m->Q, freqs, m->nstate );
    
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    //update_eigen_system( m );
    m->need_update = false;
}

void _rev_relative_update_Q( SubstitutionModel *m ){
	const double* freqs = m->simplex->get_values(m->simplex);
    int index = 0;
	const unsigned* model = m->model->values;
    for ( int i = 0; i < m->nstate; i++ ) {
        for ( int j = i+1; j < m->nstate; j++, index++ ) {
            if( i == m->nstate-2 && j == m->nstate-1 ){
                m->Q[i][j] = m->Q[j][i] = 1;
            }
            else {
                m->Q[i][j] = m->Q[j][i] = Parameters_value(m->rates, model[index]);
            }
        }
    }
    
    double temp;
    for ( int i = 0; i < m->nstate; i++ )  {
        for ( int j = i + 1; j < m->nstate; j++ ) {
            temp = m->Q[i][j];
            m->Q[i][j] = temp * freqs[j];
            m->Q[j][i] = temp * freqs[i];
        }
    }
    
    update_eigen_system( m );
    m->need_update = false;
}

/******************************************************************************************************************************************************/
#pragma mark -
#pragma mark Non Reversible


static void _genernal_update_Q( SubstitutionModel *m );


void _genernal_update_Q( SubstitutionModel *m ){
    int index = 0;
	const unsigned* model = m->model->values;
    for ( int i = 0; i < m->nstate; i++ )  {
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++])/* * freqs[j]*/;
        }
        j++;
        for ( ; j < m->nstate; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, model[index++])/* * freqs[j]*/;
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
    //normalize_Q( m->Q, freqs, m->nstate );
    
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    //update_eigen_system( m );
    m->need_update = false;
}

