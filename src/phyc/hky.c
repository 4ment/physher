/*
 *  hky.c
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

#include "hky.h"

#include "nucsubst.h"
#include "matrix.h"

/***************
	HKY
 Q matrix
 
    A C G T
 A  * . k .
 C  . * . k
 G  k . * .
 T  . k . *
 
 freqs  A C G T
 0 1 2 3
 
 rates a b c d e
 0 1 2 3 4
 1 k 1 1 k
 ***************/

static void hky_update_Q( SubstitutionModel *model );
static void _hky_p_t( SubstitutionModel *m, const double t, double *P );
static void _hky_p_t_transpose( SubstitutionModel *m, const double t, double *P );

static void _hky_dp_dt( SubstitutionModel *m, const double t, double *P );
static void _hky_dp_dt_transpose( SubstitutionModel *m, const double t, double *P );

static void _hky_d2p_dt2( SubstitutionModel *m, const double t, double *P );
static void _hky_d2p_dt2_transpose( SubstitutionModel *m, const double t, double *P );

static void _foo_update( SubstitutionModel *m ){} // does not nothing

SubstitutionModel * new_HKY(){
    double freqs[4] = {0.25, 0.25, 0.25, 0.25};
    double kappa = 3;
    return new_HKY_with_values(freqs, kappa);
}

SubstitutionModel * new_HKY_with_values( const double *freqs, const double kappa ){
	
	Simplex* sfreqs = new_Simplex(4);
	sfreqs->set_values(sfreqs, freqs);
	
	SubstitutionModel *m = create_nucleotide_model("HKY", HKY, sfreqs);
	
	m->update_Q = _foo_update;
	m->p_t = _hky_p_t;
	m->p_t_transpose = _hky_p_t_transpose;
	
	m->dp_dt = _hky_dp_dt;
	m->dp_dt_transpose = _hky_dp_dt_transpose;
	
	m->rates = new_Parameters( 1 );
	Parameters_move(m->rates, new_Parameter_with_postfix("hky.kappa", "model", kappa, new_Constraint(0.0001, 100) ) );

	return m;
}

SubstitutionModel * new_HKY_with_parameters( Simplex *freqs, const Parameters* kappa ){
	
	SubstitutionModel *m = create_nucleotide_model("HKY", HKY, freqs);
	
	m->update_Q = _foo_update;
	m->p_t = _hky_p_t;
	m->p_t_transpose = _hky_p_t_transpose;
	
	m->dp_dt = _hky_dp_dt;
	m->dp_dt_transpose = _hky_dp_dt_transpose;
	
	m->rates = new_Parameters( 1 );
	Parameters_add(m->rates, Parameters_at(kappa, 0) );
	
	return m;
}




void hky_update_Q( SubstitutionModel *model ){
	
    model->Q[1][3] = model->Q[3][1] = model->Q[0][2] = model->Q[2][0] = Parameters_value(model->rates, 0); // kappa
    model->Q[0][1] = model->Q[1][0] = model->Q[0][3] = model->Q[3][0] = model->Q[1][2] = model->Q[2][1] = model->Q[2][3] = model->Q[3][2] = 1.;
    
	const double* freqs = model->get_frequencies(model);
    double temp;
    for ( int i = 0; i < 4; i++ )  {
        for ( int j = i + 1; j < 4; j++ ) {
            temp = model->Q[i][j];
            model->Q[i][j] = temp * freqs[j];
            model->Q[j][i] = temp * freqs[i];
        }
    }
    
    update_eigen_system( model );
    model->need_update = false;
}

void hky_update_Q2( SubstitutionModel *m ){
    double kappa = Parameters_value(m->rates, 0);
	const double* freqs = m->get_frequencies(m);
	
    m->Q[0][0] =  freqs[1] + kappa*freqs[2]+ freqs[3];
    m->Q[1][1] =  freqs[0] + freqs[2]+ kappa*freqs[3];
    m->Q[2][2] =  kappa*freqs[0] + freqs[1]+ kappa*freqs[3];
    m->Q[3][3] =  freqs[0] + kappa*freqs[1]+ freqs[2];
    
    m->Q[1][0] = m->Q[3][0] = freqs[0];
    m->Q[0][1] = m->Q[1][2] = freqs[1];
    m->Q[1][2] = m->Q[3][2] = freqs[2];
    m->Q[0][3] = m->Q[1][3] = freqs[3];
    
    m->Q[2][0] = kappa * freqs[0];
    m->Q[3][1] = kappa * freqs[1];
    m->Q[0][2] = kappa * freqs[2];
    m->Q[1][3] = kappa * freqs[3];
    
    double r = 1. / (2. * (freqs[0] * freqs[1] + freqs[1] * freqs[2] + freqs[0] * freqs[3] + freqs[2] * freqs[3] + kappa * (freqs[1] * freqs[3] + freqs[0] * freqs[2])));
    
    for ( int i = 0; i < 4; i++ )  {
        for ( int j = i + 1; j < 4; j++ ) {
            m->Q[i][j] *= r;
        }
    }
    m->Q[0][0] *= freqs[0];
    m->Q[0][1] *= freqs[1];
    m->Q[0][2] *= freqs[2];
    m->Q[0][3] *= freqs[3];
    
    m->Q[1][0] *= freqs[0];
    m->Q[1][1] *= freqs[1];
    m->Q[1][2] *= freqs[2];
    m->Q[1][3] *= freqs[3];
    
    m->Q[2][0] *= freqs[0];
    m->Q[2][1] *= freqs[1];
    m->Q[2][2] *= freqs[2];
    m->Q[2][3] *= freqs[3];
    
    m->Q[3][0] *= freqs[0];
    m->Q[3][1] *= freqs[1];
    m->Q[3][2] *= freqs[2];
    m->Q[3][3] *= freqs[3];
    
    //update_eigen_system( m );
    m->need_update = false;
}

// row major
void _hky_p_t( SubstitutionModel *m, const double t, double *P ){
    
    //	if( m->need_update ){
    //		m->update_Q(m);
    //	}
    const double* freqs = m->get_frequencies(m);
    
    double R  = freqs[0] + freqs[2];
    double Y  = freqs[3] + freqs[1];
    double k1 = Parameters_value(m->rates, 0) * Y + R;
    double k2 = Parameters_value(m->rates, 0) * R + Y;
    
    double r = 1. / (2. * (freqs[0] * freqs[1] + freqs[1] * freqs[2] + freqs[0] * freqs[3] + freqs[2] * freqs[3] + Parameters_value(m->rates, 0) * (freqs[1] * freqs[3] + freqs[0] * freqs[2])));
    
    double exp1  = exp(-t*r);
    double exp22 = exp(-k2 * t * r);
    double exp21 = exp(-k1 * t * r);
    
    //A
    P[0]  = freqs[0] * (1. + (Y/R) * exp1) + (freqs[2]/R) * exp22; //A
    P[1]  = freqs[1] * (1. -               exp1);                        //C
    P[2]  = freqs[2] * (1. + (Y/R) * exp1) - (freqs[2]/R) * exp22; //G
    P[3]  = freqs[3] * (1. -               exp1);                        //T, U
    
    //C
    P[4]  = freqs[0] * (1. -               exp1);                        //A
    P[5]  = freqs[1] * (1. + (R/Y) * exp1) + (freqs[3]/Y) * exp21; //C
    P[6]  = freqs[2] * (1. -               exp1);                        //G
    P[7]  = freqs[3] * (1. + (R/Y) * exp1) - (freqs[3]/Y) * exp21; //T, U
    
    //G
    P[8]  = freqs[0] * (1. + (Y/R) * exp1) - (freqs[0]/R) * exp22; //A
    P[9]  = freqs[1] * (1. -               exp1);                        //C
    P[10] = freqs[2] * (1. + (Y/R) * exp1) + (freqs[0]/R) * exp22; //G
    P[11] = freqs[3] * (1. -               exp1);                        //T, U
    
    //T, U
    P[12] = freqs[0] * (1. -               exp1);                        //A
    P[13] = freqs[1] * (1. + (R/Y) * exp1) - (freqs[1]/Y) * exp21; //C
    P[14] = freqs[2] * (1. -               exp1);                        //G
    P[15] = freqs[3] * (1. + (R/Y) * exp1) + (freqs[1]/Y) * exp21; //T, U
    
}

// column major
void _hky_p_t_transpose( SubstitutionModel *m, const double t, double *P ){
    
    //	if( m->need_update ){
    //		m->update_Q(m);
    //	}
    const double* freqs = m->get_frequencies(m);
    
    double R  = freqs[0] + freqs[2];
    double Y  = freqs[3] + freqs[1];
    double k1 = Parameters_value(m->rates, 0) * Y + R;
    double k2 = Parameters_value(m->rates, 0) * R + Y;
    
    double r = 1. / (2. * (freqs[0] * freqs[1] + freqs[1] * freqs[2] + freqs[0] * freqs[3] + freqs[2] * freqs[3] + Parameters_value(m->rates, 0) * (freqs[1] * freqs[3] + freqs[0] * freqs[2])));
    
    double exp1  = exp(-t*r);
    double exp22 = exp(-k2 * t * r);
    double exp21 = exp(-k1 * t * r);
    
    //A
    P[0]  = freqs[0] * (1. + (Y/R) * exp1) + (freqs[2]/R) * exp22; //A
    P[4]  = freqs[1] * (1. -               exp1);                        //C
    P[8]  = freqs[2] * (1. + (Y/R) * exp1) - (freqs[2]/R) * exp22; //G
    P[12] = freqs[3] * (1. -               exp1);                        //T, U
    
    //C
    P[1]  = freqs[0] * (1. -               exp1);                        //A
    P[5]  = freqs[1] * (1. + (R/Y) * exp1) + (freqs[3]/Y) * exp21; //C
    P[9]  = freqs[2] * (1. -               exp1);                        //G
    P[13] = freqs[3] * (1. + (R/Y) * exp1) - (freqs[3]/Y) * exp21; //T, U
    
    //G
    P[2]  = freqs[0] * (1. + (Y/R) * exp1) - (freqs[0]/R) * exp22; //A
    P[6]  = freqs[1] * (1. -               exp1);                        //C
    P[10] = freqs[2] * (1. + (Y/R) * exp1) + (freqs[0]/R) * exp22; //G
    P[14] = freqs[3] * (1. -               exp1);                        //T, U
    
    //T, U
    P[3]  = freqs[0] * (1. -               exp1);                        //A
    P[7]  = freqs[1] * (1. + (R/Y) * exp1) - (freqs[1]/Y) * exp21; //C
    P[11] = freqs[2] * (1. -               exp1);                        //G
    P[15] = freqs[3] * (1. + (R/Y) * exp1) + (freqs[1]/Y) * exp21; //T, U
    
}


void _hky_dp_dt( SubstitutionModel *m, const double t, double *P ) {
    
    const double* freqs = m->get_frequencies(m);
    
    double R  = freqs[0] + freqs[2];
    double Y  = freqs[3] + freqs[1];
    double k1 = Parameters_value(m->rates, 0) * Y + R;
    double k2 = Parameters_value(m->rates, 0) * R + Y;
    
    double r = 1. / (2. * (freqs[0] * freqs[1] + freqs[1] * freqs[2] + freqs[0] * freqs[3] + freqs[2] * freqs[3] + Parameters_value(m->rates, 0) * (freqs[1] * freqs[3] + freqs[0] * freqs[2])));
    
    double exp1  = exp(-t*r);
    double exp22 = exp(-k2 * t * r);
    double exp21 = exp(-k1 * t * r);
    
    //A
    P[0]  = r * (freqs[0] * -(Y/R) * exp1 - (freqs[2]/R) * k2 * exp22); //A
    P[1]  = r * (freqs[1] *                exp1);                              //C
    P[2]  = r * (freqs[2] * -(Y/R) * exp1 + (freqs[2]/R) * k2 * exp22); //G
    P[3]  = r * (freqs[3] *                exp1);                              //T, U
    
    //C
    P[4]  = r * (freqs[0] *                exp1);                              //A
    P[5]  = r * (freqs[1] * -(R/Y) * exp1 - (freqs[3]/Y) * k1 * exp21); //C
    P[6]  = r * (freqs[2] *                exp1);                              //G
    P[7]  = r * (freqs[3] * -(R/Y) * exp1 + (freqs[3]/Y) * k1 * exp21); //T, U
    
    //G
    P[8]  = r * (freqs[0] * -(Y/R) * exp1 + (freqs[0]/R) * k2 * exp22); //A
    P[9]  = r * (freqs[1] *                exp1);                              //C
    P[10] = r * (freqs[2] * -(Y/R) * exp1 - (freqs[0]/R) * k2 * exp22); //G
    P[11] = r * (freqs[3] *                exp1);                              //T, U
    
    //T, U
    P[12] = r * (freqs[0] *                exp1);                              //A
    P[13] = r * (freqs[1] * -(R/Y) * exp1 + (freqs[1]/Y) * k1 * exp21); //C
    P[14] = r * (freqs[2] *                exp1);                              //G
    P[15] = r * (freqs[3] * -(R/Y) * exp1 - (freqs[1]/Y) * k1 * exp21); //T, U
    
    
}

void _hky_dp_dt_transpose( SubstitutionModel *m, const double t, double *P ) {
    //fprintf(stderr, "dp_dt_T %d\n", counter++);
    const double* freqs = m->get_frequencies(m);
    
    double R  = freqs[0] + freqs[2];
    double Y  = freqs[3] + freqs[1];
    double k1 = Parameters_value(m->rates, 0) * Y + R;
    double k2 = Parameters_value(m->rates, 0) * R + Y;
    
    double r = 1. / (2. * (freqs[0] * freqs[1] + freqs[1] * freqs[2] + freqs[0] * freqs[3] + freqs[2] * freqs[3] + Parameters_value(m->rates, 0) * (freqs[1] * freqs[3] + freqs[0] * freqs[2])));
    
    double exp1  = exp(-t*r);
    double exp22 = exp(-k2 * t * r);
    double exp21 = exp(-k1 * t * r);
    
    //A
    P[0]  = r * (freqs[0] * -(Y/R) * exp1 - (freqs[2]/R) * k2 * exp22); //A
    P[4]  = r * (freqs[1] *                exp1);                              //C
    P[8]  = r * (freqs[2] * -(Y/R) * exp1 + (freqs[2]/R) * k2 * exp22); //G
    P[12]  = r * (freqs[3] *                exp1);                              //T, U
    
    //C
    P[1]  = r * (freqs[0] *                exp1);                              //A
    P[5]  = r * (freqs[1] * -(R/Y) * exp1 - (freqs[3]/Y) * k1 * exp21); //C
    P[9]  = r * (freqs[2] *                exp1);                              //G
    P[13] = r * (freqs[3] * -(R/Y) * exp1 + (freqs[3]/Y) * k1 * exp21); //T, U
    
    //G
    P[2]  = r * (freqs[0] * -(Y/R) * exp1 + (freqs[0]/R) * k2 * exp22); //A
    P[6]  = r * (freqs[1] *                exp1);                              //C
    P[10] = r * (freqs[2] * -(Y/R) * exp1 - (freqs[0]/R) * k2 * exp22); //G
    P[14] = r * (freqs[3] *                exp1);                              //T, U
    
    //T, U
    P[3]  = r * (freqs[0] *                exp1);                              //A
    P[7]  = r * (freqs[1] * -(R/Y) * exp1 + (freqs[1]/Y) * k1 * exp21); //C
    P[11] = r * (freqs[2] *                exp1);                              //G
    P[15] = r * (freqs[3] * -(R/Y) * exp1 - (freqs[1]/Y) * k1 * exp21); //T, U
    //print_dvector(P, 16);exit(1);
}

void _hky_d2p_dt2( SubstitutionModel *m, const double t, double *P ) {
    
    const double* freqs = m->get_frequencies(m);
    
    double R  = freqs[0] + freqs[2];
    double Y  = freqs[3] + freqs[1];
    double k1 = Parameters_value(m->rates, 0) * Y + R;
    double k2 = Parameters_value(m->rates, 0) * R + Y;
    double k12 = k1*k1;
    double k22 = k2*k2;
    
    double r = 1. / (2. * (freqs[0] * freqs[1] + freqs[1] * freqs[2] + freqs[0] * freqs[3] + freqs[2] * freqs[3] + Parameters_value(m->rates, 0) * (freqs[1] * freqs[3] + freqs[0] * freqs[2])));
    double r2 = r * r;
    
    double l = r * t;
    double exp1  = exp(-l);
    double exp22 = exp(-k2 * l);
    double exp21 = exp(-k1 * l);
    
    //A
    P[0] = r2 * (freqs[0] * (Y/R) * exp1 + (freqs[2]/R) * k22 * exp22); //A
    P[1] = r2 * (freqs[1] *             - exp1);                               //C
    P[2] = r2 * (freqs[2] * (Y/R) * exp1 - (freqs[2]/R) * k22 * exp22); //G
    P[3] = r2 * (freqs[3] *             - exp1);                               //T, U
    
    //C
    P[4] = r2 * (freqs[0] *             - exp1);                               //A
    P[5] = r2 * (freqs[1] * (R/Y) * exp1 + (freqs[3]/Y) * k12 * exp21); //C
    P[6] = r2 * (freqs[2] *             - exp1);                               //G
    P[7] = r2 * (freqs[3] * (R/Y) * exp1 - (freqs[3]/Y) * k12 * exp21); //T, U
    
    //G
    P[8] = r2 * (freqs[0] * (Y/R) * exp1 - (freqs[0]/R) * k22 * exp22); //A
    P[9] = r2 * (freqs[1] *             - exp1);                               //C
    P[10] = r2 * (freqs[2] * (Y/R) * exp1 + (freqs[0]/R) * k22 * exp22); //G
    P[11] = r2 * (freqs[3] *             - exp1);                               //T, U
    
    //T, U
    P[12] = r2 * (freqs[0] *             - exp1);                               //A
    P[13] = r2 * (freqs[1] * (R/Y) * exp1 - (freqs[1]/Y) * k12 * exp21); //C
    P[14] = r2 * (freqs[2] *             - exp1);                               //G
    P[15] = r2 * (freqs[3] * (R/Y) * exp1 + (freqs[1]/Y) * k12 * exp21); //T, U
    
    
}

void _hky_d2p_dt2_transpose( SubstitutionModel *m, const double t, double *P ) {
    
    const double* freqs = m->get_frequencies(m);
    
    double R  = freqs[0] + freqs[2];
    double Y  = freqs[3] + freqs[1];
    double k1 = Parameters_value(m->rates, 0) * Y + R;
    double k2 = Parameters_value(m->rates, 0) * R + Y;
    double k12 = k1*k1;
    double k22 = k2*k2;
    
    double r = 1. / (2. * (freqs[0] * freqs[1] + freqs[1] * freqs[2] + freqs[0] * freqs[3] + freqs[2] * freqs[3] + Parameters_value(m->rates, 0) * (freqs[1] * freqs[3] + freqs[0] * freqs[2])));
    double r2 = r * r;
    
    double l = r * t;
    double exp1  = exp(-l);
    double exp22 = exp(-k2 * l);
    double exp21 = exp(-k1 * l);
    
    //A
    P[0] = r2 * (freqs[0] * (Y/R) * exp1 + (freqs[2]/R) * k22 * exp22); //A
    P[4] = r2 * (freqs[1] *             - exp1);                               //C
    P[8] = r2 * (freqs[2] * (Y/R) * exp1 - (freqs[2]/R) * k22 * exp22); //G
    P[12] = r2 * (freqs[3] *             - exp1);                               //T, U
    
    //C
    P[1] = r2 * (freqs[0] *             - exp1);                               //A
    P[5] = r2 * (freqs[1] * (R/Y) * exp1 + (freqs[3]/Y) * k12 * exp21); //C
    P[9] = r2 * (freqs[2] *             - exp1);                               //G
    P[13] = r2 * (freqs[3] * (R/Y) * exp1 - (freqs[3]/Y) * k12 * exp21); //T, U
    
    //G
    P[2] = r2 * (freqs[0] * (Y/R) * exp1 - (freqs[0]/R) * k22 * exp22); //A
    P[6] = r2 * (freqs[1] *             - exp1);                               //C
    P[10] = r2 * (freqs[2] * (Y/R) * exp1 + (freqs[0]/R) * k22 * exp22); //G
    P[14] = r2 * (freqs[3] *             - exp1);                               //T, U
    
    //T, U
    P[3] = r2 * (freqs[0] *             - exp1);                               //A
    P[7] = r2 * (freqs[1] * (R/Y) * exp1 - (freqs[1]/Y) * k12 * exp21); //C
    P[11] = r2 * (freqs[2] *             - exp1);                               //G
    P[15] = r2 * (freqs[3] * (R/Y) * exp1 + (freqs[1]/Y) * k12 * exp21); //T, U
}
