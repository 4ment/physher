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


bool isSymmetric(const unsigned *matrix, unsigned dim){
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
    int index = 0;
    for ( int i = 0; i < m->nstate; i++ )  {
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, m->model[index++]) * m->_freqs[j];
        }
        j++;
        for ( ; j < m->nstate; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, m->model[index++]) * m->_freqs[j];
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
    if(m->normalize)normalize_Q( m->Q, m->_freqs, m->nstate );
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    
    m->need_update = false;
}

void _nonreversible_update_Q2( SubstitutionModel *m ){
    int index = 0;
    for ( int i = 0; i < m->nstate; i++ )  {
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, m->model[index++]) * m->_freqs[j];
        }
        j++;
        for ( ; j < m->nstate; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, m->model[index++]) * m->_freqs[j];
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
    if(m->normalize)normalize_Q( m->Q, m->_freqs, m->nstate );
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    
    m->need_update = false;
}


void _reversible_update_Q( SubstitutionModel *m ){
    double temp;
    for ( int i = 0; i < m->nstate; i++ )  {
        for ( int j = i + 1; j < m->nstate; j++ ) {
            temp = Parameters_value(m->rates, m->model[i*m->nstate+j]);
            m->Q[i][j] = temp * m->_freqs[j];
            m->Q[j][i] = temp * m->_freqs[i];
        }
    }
    make_zero_rows( m->Q, m->nstate);
    if(m->normalize)normalize_Q( m->Q, m->_freqs, m->nstate );
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    
    m->need_update = false;
}

SubstitutionModel * new_GeneralModel( const unsigned *model, unsigned dim ){
    SubstitutionModel *m = new_GeneralModel2(model,dim,-1,true);
    return m;
}

SubstitutionModel * new_GeneralModel2( const unsigned *model, unsigned dim, int relativeTo, bool normalize ){
    bool sym = isSymmetric(model, dim);
    SubstitutionModel *m = NULL;
    if(sym){
        m = create_general_model("GENERAL", REVERSIBLE, dim);
        m->update_Q = _reversible_update_Q;
    }
    else{
        m = create_general_model("UREVGENERAL", NONREVERSIBLE, dim);
        m->update_Q = _nonreversible_update_Q;
        m->reversible = false;
    }
    m->normalize = normalize;
    
    m->_freqs = dvector(dim);
    
    for ( int i = 0; i < dim; i++ ) {
        m->_freqs[i] = 1.0/dim;
    }
    
    int len = dim*dim;
    m->model = clone_uivector(model, len);
    
    int max = uimax_vector(model, len);
    max++;
    
    
    StringBuffer *buffer = new_StringBuffer(10);
    
    m->rates = new_Parameters( max);
    for ( int i = 0; i < max; i++ ) {
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "q.%d", i);
        Parameters_move(m->rates, new_Parameter_with_postfix(buffer->c, "model", 1.0, new_Constraint(0.00001, 1000) ) );
    }
    
    free_StringBuffer(buffer);
    
    return m;
}

#pragma mark -
#pragma mark Reversible


static void _rev_update_Q( SubstitutionModel *m );
static void _rev_relative_update_Q( SubstitutionModel *m );

SubstitutionModel * new_ReversibleModel( int dim, const int *model ){
    SubstitutionModel *m = NULL;
    
    m = create_general_model("reversible", REVERSIBLE, dim);
    m->_freqs = dvector(dim);
    //default model has fixed frequencies
    for ( int i = 0; i < dim; i++ ) {
        m->_freqs[i] = 1.0/dim;
    }
    
    int len = dim*(dim-1)/2;
    m->model = uivector(len);
    
    int max = 0;
    memcpy(m->model, model, len*sizeof(int));
    max = imax_vector(model, len);
    max++;
    
    m->update_Q = _rev_update_Q;
    
    StringBuffer *buffer = new_StringBuffer(10);
    m->rates = new_Parameters( max);
    for ( int i = 0; i < max; i++ ) {
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "q.%d", i);
        Parameters_move(m->rates, new_Parameter_with_postfix(buffer->c, "model", 1.0, new_Constraint(0.00001, 1000) ) );
    }
    
    free_StringBuffer(buffer);
    
    return m;
}

void _rev_update_Q( SubstitutionModel *m ){
    int index = 0;
    double temp;
    for ( int i = 0; i < m->nstate; i++ )  {
        for ( int j = i + 1; j < m->nstate; j++ ) {
            temp = Parameters_value(m->rates, m->model[index++]);
            m->Q[i][j] = temp;// * m->_freqs[j];
            m->Q[j][i] = temp;// * m->_freqs[i];
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
    //normalize_Q( m->Q, m->_freqs, m->nstate );
    
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    //update_eigen_system( m );
    m->need_update = false;
}

void _rev_relative_update_Q( SubstitutionModel *m ){
    int index = 0;
    for ( int i = 0; i < m->nstate; i++ ) {
        for ( int j = i+1; j < m->nstate; j++, index++ ) {
            if( i == m->nstate-2 && j == m->nstate-1 ){
                m->Q[i][j] = m->Q[j][i] = 1;
            }
            else {
                m->Q[i][j] = m->Q[j][i] = Parameters_value(m->rates, m->model[index]);
            }
        }
    }
    
    double temp;
    for ( int i = 0; i < m->nstate; i++ )  {
        for ( int j = i + 1; j < m->nstate; j++ ) {
            temp = m->Q[i][j];
            m->Q[i][j] = temp * m->_freqs[j];
            m->Q[j][i] = temp * m->_freqs[i];
        }
    }
    
    update_eigen_system( m );
    m->need_update = false;
}

/******************************************************************************************************************************************************/
#pragma mark -
#pragma mark Non Reversible


static void _genernal_update_Q( SubstitutionModel *m );

SubstitutionModel * new_NonReversibleModel( int dim, const int *model ){
    SubstitutionModel *m = NULL;
    
    m = create_general_model("non.reversible", NONREVERSIBLE, dim);
    m->_freqs = dvector(dim);
    //default model has fixed equal frequencies
    for ( int i = 0; i < dim; i++ ) {
        m->_freqs[i] = 1.0/dim;
    }
    
    int len = dim*dim-dim;
    m->model = uivector(len);
    memcpy(m->model, model, len*sizeof(int));
    int max = 0;
    max = imax_vector(model, len);
    max++;
    
    m->update_Q = _genernal_update_Q;
    
    StringBuffer *buffer = new_StringBuffer(10);
    m->rates = new_Parameters( max);
    for ( int i = 0; i < max; i++ ) {
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "q.%d", i);
        Parameters_move(m->rates, new_Parameter_with_postfix(buffer->c, "model", 1.0, new_Constraint(0.00001, 1000) ) );
    }
    
    free_StringBuffer(buffer);
    
    return m;
}

void _genernal_update_Q( SubstitutionModel *m ){
    int index = 0;
    for ( int i = 0; i < m->nstate; i++ )  {
        int j = 0;
        for ( ; j < i; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, m->model[index++])/* * m->_freqs[j]*/;
        }
        j++;
        for ( ; j < m->nstate; j++ ) {
            m->Q[i][j] = Parameters_value(m->rates, m->model[index++])/* * m->_freqs[j]*/;
        }
    }
    
    make_zero_rows( m->Q, m->nstate);
    //normalize_Q( m->Q, m->_freqs, m->nstate );
    
    EigenDecomposition_decompose(m->Q, m->eigendcmp);
    //update_eigen_system( m );
    m->need_update = false;
}

