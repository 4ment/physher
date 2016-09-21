/*
 *  f81.h
 *  physher
 *
 *  Created by Mathieu Fourment on 19/09/2016.
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

#include "f81.h"

#include "matrix.h"

static void f81_p_t( SubstitutionModel *m, const double t, double *P );
static void f81_p_t_transpose( SubstitutionModel *m, const double t, double *P );
static double f81_pij_t( SubstitutionModel *m, const int i, const int j, const double t );

static void f81_dp_dt( SubstitutionModel *m, const double t, double *P );
static void f81_dp_dt_transpose( SubstitutionModel *m, const double t, double *P );
static void f81_d2p_dt2( SubstitutionModel *m, const double t, double *P );
static void f81_d2p_dt2_transpose( SubstitutionModel *m, const double t, double *P );

SubstitutionModel * new_F81(){
    
    SubstitutionModel *m = create_substitution_model("F81", JC69, DATA_TYPE_NUCLEOTIDE);
    m->nstate = 4;
    m->_freqs = dvector(4);
    for ( int i = 0; i < 4; i++ ) m->_freqs[i] = 0.25;
    
    m->freqs = new_Parameters( 3 );
    m->update_frequencies = nucleotide_update_freqs;
    double aux1 = m->_freqs[1] / (  1 - m->_freqs[0]);
    double aux2 = m->_freqs[2] / ( (1 - m->_freqs[0]) * (1 - aux1) );
    Parameters_add(m->freqs, new_Parameter_with_postfix("hky.piA", "model", m->_freqs[0], new_Constraint(0.001, 0.999) ) );
    Parameters_add(m->freqs, new_Parameter_with_postfix("hky.piC", "model", aux1,     new_Constraint(0.001, 0.999) ) );
    Parameters_add(m->freqs, new_Parameter_with_postfix("hky.piT", "model", aux2,     new_Constraint(0.001, 0.999) ) );
    
    m->pij_t = f81_pij_t;
    m->p_t = f81_p_t;
    m->p_t_transpose = f81_p_t_transpose;
    m->dp_dt = f81_dp_dt;
    m->dp_dt_transpose = f81_dp_dt_transpose;
    m->d2p_d2t = f81_d2p_dt2;
    m->d2p_d2t_transpose = f81_d2p_dt2_transpose;
    return m;
}

double f81_pij_t( SubstitutionModel *m, const int i, const int j, const double t ){
    double temp =  exp(-t/(1.0 - m->_freqs[0]*m->_freqs[0] - m->_freqs[1]*m->_freqs[1] - m->_freqs[2]*m->_freqs[2] - m->_freqs[3]*m->_freqs[3]));
    if(i == j) return temp + m->_freqs[j]*(1.0-temp);
    else       return m->_freqs[j]*(1.0-temp);
}

void f81_p_t_rolled( SubstitutionModel *m, const double t, double *P ){
    double temp =  exp(-t/(1.0 - m->_freqs[0]*m->_freqs[0] - m->_freqs[1]*m->_freqs[1] - m->_freqs[2]*m->_freqs[2] - m->_freqs[3]*m->_freqs[3]));
    int i,j;
    int k = 0;
    for ( i = 0; i < 4; i++) {
        for ( j = 0; j < 4; j++){
            if( i == j ) P[k] = temp + m->_freqs[j]*(1.0-temp);
            else P[k] = m->_freqs[j]*(1.0-temp);
            k++;
        }
    }
}

void f81_p_t( SubstitutionModel *m, const double t, double *P ){
    double temp =  exp(-t/(1.0 - m->_freqs[0]*m->_freqs[0] - m->_freqs[1]*m->_freqs[1] - m->_freqs[2]*m->_freqs[2] - m->_freqs[3]*m->_freqs[3]));
    double temp2 = 1.0 - temp;
    
    double atemp2 = m->_freqs[0] * temp2;
    double ctemp2 = m->_freqs[1] * temp2;
    double ttemp2 = m->_freqs[2] * temp2;
    double gtemp2 = m->_freqs[3] * temp2;
    
    //A
    P[0]  = atemp2 + temp; //A
    P[1]  = ctemp2;        //C
    P[2]  = ttemp2;        //G
    P[3]  = gtemp2;        //T, U
    
    //C
    P[4]  = atemp2;        //A
    P[5]  = ctemp2 + temp; //C
    P[6]  = ttemp2;        //G
    P[7]  = gtemp2;        //T, U
    
    //G
    P[8]  = atemp2;        //A
    P[9]  = ctemp2;        //C
    P[10] = ttemp2 + temp; //G
    P[11] = gtemp2;        //T, U
    
    //T, U
    P[12] = atemp2;        //A
    P[13] = ctemp2;        //C
    P[14] = ttemp2;        //G
    P[15] = gtemp2 + temp; //T, U
}

void f81_p_t_transpose( SubstitutionModel *m, const double t, double *P ){
    double temp =  exp(-t/(1.0 - m->_freqs[0]*m->_freqs[0] - m->_freqs[1]*m->_freqs[1] - m->_freqs[2]*m->_freqs[2] - m->_freqs[3]*m->_freqs[3]));
    double temp2 = 1.0 - temp;
    
    double atemp2 = m->_freqs[0] * temp2;
    double ctemp2 = m->_freqs[1] * temp2;
    double ttemp2 = m->_freqs[2] * temp2;
    double gtemp2 = m->_freqs[3] * temp2;
    
    //A
    P[0]  = atemp2 + temp; //A
    P[4]  = ctemp2;        //C
    P[8]  = ttemp2;        //G
    P[12]  = gtemp2;        //T, U
    
    //C
    P[1]  = atemp2;        //A
    P[5]  = ctemp2 + temp; //C
    P[9]  = ttemp2;        //G
    P[13]  = gtemp2;        //T, U
    
    //G
    P[2]  = atemp2;        //A
    P[6]  = ctemp2;        //C
    P[10] = ttemp2 + temp; //G
    P[14] = gtemp2;        //T, U
    
    //T, U
    P[3] = atemp2;        //A
    P[7] = ctemp2;        //C
    P[11] = ttemp2;        //G
    P[15] = gtemp2 + temp; //T, U
}

void f81_dp_dt( SubstitutionModel *m, const double t, double *P ){
    double beta = 1.0/(1.0 - m->_freqs[0]*m->_freqs[0] - m->_freqs[1]*m->_freqs[1] - m->_freqs[2]*m->_freqs[2] - m->_freqs[3]*m->_freqs[3]);
    double betatemp = beta*exp(-beta*t);
    
    double atemp = m->_freqs[0] * betatemp;
    double ctemp = m->_freqs[1] * betatemp;
    double ttemp = m->_freqs[2] * betatemp;
    double gtemp = m->_freqs[3] * betatemp;
    
    //A
    P[0]  = betatemp - atemp; //A
    P[1]  = ctemp;            //C
    P[2]  = ttemp;            //G
    P[3]  = gtemp;            //T, U
    
    //C
    P[4]  = atemp;            //A
    P[5]  = betatemp - ctemp; //C
    P[6]  = ttemp;            //G
    P[7]  = gtemp;            //T, U
    
    //G
    P[8]  = atemp;            //A
    P[9]  = ctemp;            //C
    P[10] = betatemp - ttemp; //G
    P[11] = gtemp;            //T, U
    
    //T, U
    P[12] = atemp;            //A
    P[13] = ctemp;            //C
    P[14] = ttemp;            //G
    P[15] = betatemp - gtemp; //T, U
}

void f81_dp_dt_transpose( SubstitutionModel *m, const double t, double *P ){
    double beta = 1.0/(1.0 - m->_freqs[0]*m->_freqs[0] - m->_freqs[1]*m->_freqs[1] - m->_freqs[2]*m->_freqs[2] - m->_freqs[3]*m->_freqs[3]);
    double betatemp = beta*exp(-beta*t);
    
    double atemp = m->_freqs[0] * betatemp;
    double ctemp = m->_freqs[1] * betatemp;
    double ttemp = m->_freqs[2] * betatemp;
    double gtemp = m->_freqs[3] * betatemp;
    
    //A
    P[0]  = betatemp - atemp; //A
    P[4]  = ctemp;            //C
    P[8]  = ttemp;            //G
    P[12] = gtemp;            //T, U
    
    //C
    P[1]  = atemp;            //A
    P[5]  = betatemp - ctemp; //C
    P[9]  = ttemp;            //G
    P[13] = gtemp;            //T, U
    
    //G
    P[2]  = atemp;            //A
    P[6]  = ctemp;            //C
    P[10] = betatemp - ttemp; //G
    P[14] = gtemp;            //T, U
    
    //T, U
    P[3]  = atemp;            //A
    P[7]  = ctemp;            //C
    P[11] = ttemp;            //G
    P[15] = betatemp - gtemp; //T, U
}

void f81_d2p_dt2( SubstitutionModel *m, const double t, double *P ){
    double beta = 1.0/(1.0 - m->_freqs[0]*m->_freqs[0] - m->_freqs[1]*m->_freqs[1] - m->_freqs[2]*m->_freqs[2] - m->_freqs[3]*m->_freqs[3]);
    double beta2temp = beta*beta*exp(-beta*t);
    
    double atemp = m->_freqs[0] * beta2temp;
    double ctemp = m->_freqs[1] * beta2temp;
    double ttemp = m->_freqs[2] * beta2temp;
    double gtemp = m->_freqs[3] * beta2temp;
    
    //A
    P[0]  = atemp - beta2temp; //A
    P[1]  = ctemp;            //C
    P[2]  = ttemp;            //G
    P[3]  = gtemp;            //T, U
    
    //C
    P[4]  = atemp;            //A
    P[5]  = ctemp - beta2temp; //C
    P[6]  = ttemp;            //G
    P[7]  = gtemp;            //T, U
    
    //G
    P[8]  = atemp;            //A
    P[9]  = ctemp;            //C
    P[10] = ttemp - beta2temp; //G
    P[11] = gtemp;            //T, U
    
    //T, U
    P[12] = atemp;            //A
    P[13] = ctemp;            //C
    P[14] = ttemp;            //G
    P[15] = gtemp - beta2temp; //T, U
}



void f81_d2p_dt2_transpose( SubstitutionModel *m, const double t, double *P ){
    double beta = 1.0/(1.0 - m->_freqs[0]*m->_freqs[0] - m->_freqs[1]*m->_freqs[1] - m->_freqs[2]*m->_freqs[2] - m->_freqs[3]*m->_freqs[3]);
    double beta2temp = beta*beta*exp(-beta*t);
    
    double atemp = m->_freqs[0] * beta2temp;
    double ctemp = m->_freqs[1] * beta2temp;
    double ttemp = m->_freqs[2] * beta2temp;
    double gtemp = m->_freqs[3] * beta2temp;
    
    //A
    P[0]  = atemp - beta2temp; //A
    P[4]  = ctemp;            //C
    P[8]  = ttemp;            //G
    P[12] = gtemp;            //T, U
    
    //C
    P[1]  = atemp;            //A
    P[5]  = ctemp - beta2temp; //C
    P[9]  = ttemp;            //G
    P[13] = gtemp;            //T, U
    
    //G
    P[2]  = atemp;            //A
    P[6]  = ctemp;            //C
    P[10] = ttemp - beta2temp; //G
    P[14] = gtemp;            //T, U
    
    //T, U
    P[3]  = atemp;            //A
    P[7]  = ctemp;            //C
    P[11] = ttemp;            //G
    P[15] = gtemp - beta2temp; //T, U
}
