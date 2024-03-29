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
 A  * 1 k 1
 C  1 * 1 k
 G  k 1 * 1
 T  1 k 1 *
 
 freqs  A C G T
 0 1 2 3
 
 rate kappa
 ***************/

static void hky_update_Q( SubstitutionModel *model );
static void _hky_p_t( SubstitutionModel *m, const double t, double *P );
static void _hky_p_t_transpose( SubstitutionModel *m, const double t, double *P );

static void _hky_dp_dt( SubstitutionModel *m, const double t, double *P );
static void _hky_dp_dt_transpose( SubstitutionModel *m, const double t, double *P );

static void _hky_d2p_dt2( SubstitutionModel *m, const double t, double *P );
static void _hky_d2p_dt2_transpose( SubstitutionModel *m, const double t, double *P );

static void _hky_dPdp(SubstitutionModel *m, int index, double* mat, double t);

SubstitutionModel * new_HKY(Simplex* freqs){
	Parameter* kappa = new_Parameter_with_postfix("hky.kappa", "model", 3, new_Constraint(0.0001, 100) );
	SubstitutionModel* model = new_HKY_with_parameters(freqs, kappa);
	free_Parameter(kappa);
	return model;
}

SubstitutionModel * new_HKY_with_values( const double *freqs, const double kappa ){
	
	Simplex* sfreqs = new_Simplex("hky", 4);
	sfreqs->set_values(sfreqs, freqs);
	
	SubstitutionModel *m = create_nucleotide_model("HKY", HKY, sfreqs);
	
	m->update_Q = hky_update_Q;
	m->p_t = _hky_p_t;
	m->p_t_transpose = _hky_p_t_transpose;
	
	m->dp_dt = _hky_dp_dt;
	m->dp_dt_transpose = _hky_dp_dt_transpose;
	
	m->dPdp = _hky_dPdp;
	m->dQ = dvector(16);
	
	m->rates = new_Parameters( 1 );
	Parameters_move(m->rates, new_Parameter_with_postfix("hky.kappa", "model", kappa, new_Constraint(0, INFINITY) ) );
	return m;
}

SubstitutionModel * new_HKY_with_parameters( Simplex *freqs, Parameter* kappa ){
	
	SubstitutionModel *m = create_nucleotide_model("HKY", HKY, freqs);
	
	m->update_Q = hky_update_Q;
	m->p_t = _hky_p_t;
	m->p_t_transpose = _hky_p_t_transpose;
	
	m->dp_dt = _hky_dp_dt;
	m->dp_dt_transpose = _hky_dp_dt_transpose;
	
	m->dPdp = _hky_dPdp;
	m->dQ = dvector(16);
	
	m->rates = new_Parameters( 1 );
	Parameters_add(m->rates, kappa);
	
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
    
    make_zero_rows( model->Q, 4);
	model->norm = normalize_Q( model->Q, model->get_frequencies(model), 4 );
}

void _hky_update_eigen(SubstitutionModel* m){
	double kappa = Parameters_value(m->rates, 0);
	
	double pi_a = m->simplex->get_value(m->simplex, 0);
	double pi_c = m->simplex->get_value(m->simplex, 1);
	double pi_g = m->simplex->get_value(m->simplex, 2);
	double pi_t = m->simplex->get_value(m->simplex, 3);
	
	double pi_r = pi_a + pi_g;
	double pi_y = pi_c + pi_t;
	double *v = m->eigendcmp->eval;
	double **evec = m->eigendcmp->evec;
	double **invevec = m->eigendcmp->Invevec;
	
	double beta = -1.0 / (2.0 * (pi_r*pi_y + kappa * (pi_a*pi_g + pi_c*pi_t)));
	v[0] = 0;
	v[1] = beta;
	v[2] = beta*(1 + pi_y*(kappa - 1));
	v[3] = beta*(1 + pi_r*(kappa - 1));
	
	for(int i = 0; i < 4; i++){
		memset(invevec[i], 0, sizeof(double)*4);
		memset(evec[i], 0, sizeof(double)*4);
	}

	invevec[0][0] = pi_a;
	invevec[0][1] = pi_c;
	invevec[0][2] = pi_g;
	invevec[0][3] = pi_t;
	
	invevec[1][0] =  pi_a*pi_y;
	invevec[1][1] = -pi_c*pi_r;
	invevec[1][2] =  pi_g*pi_y;
	invevec[1][3] = -pi_t*pi_r;
	
	invevec[2][1] =  1;
	invevec[2][3] = -1;

	invevec[3][0] =  1;
	invevec[3][2] = -1;

	evec[0][0] =  1;
	evec[1][0] =  1;
	evec[2][0] =  1;
	evec[3][0] =  1;

	evec[0][1]  =  1./pi_r;
	evec[1][1]  = -1./pi_y;
	evec[2][1]  =  1./pi_r;
	evec[3][1]  = -1./pi_y;

	evec[1][2]  =  pi_t/pi_y;
	evec[3][2]  = -pi_c/pi_y;

	evec[0][3] =  pi_g/pi_r;
	evec[2][3] = -pi_a/pi_r;
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

static void _hky_dQdp(SubstitutionModel *m, size_t index){
	double *dQ = m->dQ;
	size_t stateCount = m->nstate;
	memset(dQ, 0.0, sizeof(double)*stateCount*stateCount);
	
	if(m->need_update){
		hky_update_Q(m);
		_hky_update_eigen(m);
		m->need_update = false;
	}

	const double* frequencies = m->simplex->get_values(m->simplex);
	double norm = m->norm;
	double dnorm = 0;
	if(index == 0){
		dQ[0]  = -frequencies[2];
		dQ[2]  =  frequencies[2];
		dQ[5]  = -frequencies[3];
		dQ[7]  =  frequencies[3];
		dQ[8]  =  frequencies[0];
		dQ[10] = -frequencies[0];
		dQ[13] =  frequencies[1];
		dQ[15] = -frequencies[1];
		dnorm = normalizing_constant_Q_flat(dQ, frequencies, 4);
	}
	else{
		dQ[2] = dQ[7] = dQ[8] = dQ[13] = Parameters_value(m->rates, 0); // kappa
		dQ[1] = dQ[3] = dQ[4] = dQ[6] = dQ[9] = dQ[11] = dQ[12] = dQ[14] = 1.;
		double dF[4];
		if(m->grad_wrt_reparam){
			m->simplex->gradient(m->simplex, index-1, dF);
		}
		else{
			memset(dF, 0.0, sizeof(double)*4);
			dF[index-1] = 1.0;
		} 
		build_Q_flat(dQ, dF, stateCount);
		
		for (size_t i = 0; i < stateCount; i++) {
			dnorm -= dQ[i*stateCount+i]*frequencies[i] + m->Q[i][i]*norm*dF[i];
		}
	}

	for(size_t i = 0; i < stateCount; i++){
		for(size_t j = 0; j < stateCount; j++){
			// (dQ . norm - norm_dQ . Q)/norm^2
			dQ[i*stateCount + j] = (dQ[i*stateCount + j] - m->Q[i][j]*dnorm)/norm;
		}
	}
}

void _hky_dPdp(SubstitutionModel *m, int index, double* mat, double t){
	if(m->need_update){
		hky_update_Q(m);
		_hky_update_eigen(m);
		m->need_update = false;
		m->dQ_need_update = true;
	}
	if(m->dQ_need_update){
		_hky_dQdp(m, index);
		m->dQ_need_update = false;
	}
	dPdp_with_dQdp(m, m->dQ, mat, t);
}

static void _hky_dPdp_old(SubstitutionModel *m, int index, double* mat, double t){
	
	double *x = m->dQ;
	
	const double kappa = Parameters_value(m->rates, 0);
	
	const double pi_a = m->simplex->get_value(m->simplex, 0);
	const double pi_c = m->simplex->get_value(m->simplex, 1);
	const double pi_g = m->simplex->get_value(m->simplex, 2);
	const double pi_t = m->simplex->get_value(m->simplex, 3);
	
	const double norm = 2. * (pi_a*pi_c + pi_c*pi_g + pi_a*pi_t + pi_g*pi_t + kappa*(pi_c*pi_t + pi_a*pi_g));
	
	if(m->need_update){
//		m->update_Q(m);
		double pi_r = pi_a + pi_g;
		double pi_y = pi_c + pi_t;
		double *v = m->eigendcmp->eval;
		double **evec = m->eigendcmp->evec;
		double **invevec = m->eigendcmp->Invevec;
		
		double beta = -1.0 / (2.0 * (pi_r*pi_y + kappa * (pi_a*pi_g + pi_c*pi_t)));
		v[0] = 0;
		v[1] = beta;
		v[2] = beta*(1 + pi_y*(kappa - 1));
		v[3] = beta*(1 + pi_r*(kappa - 1));
		
		for(int i = 0; i < 4; i++){
			memset(invevec[i], 0, sizeof(double)*4);
			memset(evec[i], 0, sizeof(double)*4);
		}

		invevec[0][0] = pi_a;
		invevec[0][1] = pi_c;
		invevec[0][2] = pi_g;
		invevec[0][3] = pi_t;
		
		invevec[1][0] =  pi_a*pi_y;
		invevec[1][1] = -pi_c*pi_r;
		invevec[1][2] =  pi_g*pi_y;
		invevec[1][3] = -pi_t*pi_r;
		
		invevec[2][1] =  1;
		invevec[2][3] = -1;

		invevec[3][0] =  1;
		invevec[3][2] = -1;

		evec[0][0] =  1;
		evec[1][0] =  1;
		evec[2][0] =  1;
		evec[3][0] =  1;

		evec[0][1]  =  1./pi_r;
		evec[1][1]  = -1./pi_y;
		evec[2][1]  =  1./pi_r;
		evec[3][1]  = -1./pi_y;

		evec[1][2]  =  pi_t/pi_y;
		evec[3][2]  = -pi_c/pi_y;

		evec[0][3] =  pi_g/pi_r;
		evec[2][3] = -pi_a/pi_r;
		
		m->need_update = false;
	}
	
	if(m->dQ_need_update){
		
		const double phi1 = Parameters_value(m->simplex->parameters, 0);
		const double phi2 = Parameters_value(m->simplex->parameters, 1);
		const double phi3 = Parameters_value(m->simplex->parameters, 2);
		
		const double norm2 = norm*norm;
		
		memset(x, 0, sizeof(double)*16);
		// kappa
		if(index == 0){
			double dnorm = 2.*(pi_c*pi_t + pi_a*pi_g);
			double ratio = dnorm/norm2;
			x[0] = -pi_g/norm - (-pi_t - kappa*pi_g - pi_c)*ratio;
			x[1] = -pi_c*ratio;
			x[2] = pi_g/norm - (kappa*pi_g)*ratio;
			x[3] = -pi_t*ratio;
			
			x[4] = -pi_a*ratio;
			x[5] = -pi_t/norm - (-kappa*pi_t - pi_g - pi_a)*ratio;
			x[6] = -pi_g*ratio;
			x[7] = pi_t/norm - (kappa*pi_t)*ratio;
			
			x[8] = pi_a/norm - (kappa*pi_a)*ratio;
			x[9] = -pi_c*ratio;
			x[10] = -pi_a/norm - (-kappa*pi_a - pi_c - pi_t)*ratio;
			x[11] = -pi_t*ratio;
			
			x[12] = -pi_a*ratio;
			x[13] = pi_c/norm - (kappa*pi_c)*ratio;
			x[14] = -pi_g*ratio;
			x[15] = -pi_c/norm - (-kappa*pi_c - pi_g - pi_a)*ratio;
		}
		// dQdphi1
		else if(index == 1){
			double dnorm = -4*(kappa*phi1 - kappa)*phi2*phi2 - 4*((phi1 - 1)*phi2*phi2 - 2*(phi1 - 1)*phi2 + phi1 - 1)*phi3*phi3 + 4*(kappa*phi1 - kappa)*phi2 + 2*(2*(kappa*phi1 - kappa)*phi2*phi2 - 2*(kappa - 2)*phi1 + (kappa - 4*phi1 + 3)*phi2 + kappa - 3)*phi3 - 4*phi1 + 2;
			double ratio = dnorm/norm2;
			x[0] = (-kappa*(phi2 - 1)*phi3 + (phi2 - 1)*phi3 + 1)/norm - ratio*(-pi_t - kappa*pi_g - pi_c);
			x[1] = (-phi2)/norm - ratio*pi_c;
			x[2] = (kappa*(phi2 - 1)*phi3)/norm - ratio*kappa*pi_g;
			x[3] = (-((phi2 - 1)*phi3 - phi2 + 1))/norm - ratio*pi_t;
			
			x[4] = 1/norm - ratio*pi_a;
			x[5] = (((phi2 - 1)*phi3 - phi2 + 1)*kappa - (phi2 - 1)*phi3 - 1)/norm - ratio*(-kappa*pi_t - pi_g - pi_a);
			x[6] = (phi2 - 1)*phi3/norm - ratio*pi_g;
			x[7] = (-((phi2 - 1)*phi3 - phi2 + 1)*kappa)/norm - ratio*kappa*pi_t;
			
			x[8] = kappa/norm - ratio*kappa*pi_a;
			x[9] = -phi2/norm - ratio*pi_c;
			x[10] = ((phi2 - 1)*phi3 - kappa + 1)/norm - ratio*(-kappa*pi_a - pi_c - pi_t);
			x[11] = (-(phi2 - 1)*phi3 + phi2 - 1)/norm - ratio*pi_t;
			
			x[12] = 1/norm - ratio*pi_a;
			x[13] = -kappa*phi2/norm - ratio*kappa*pi_c;
			x[14] = (phi2 - 1)*phi3/norm - ratio*pi_g;
			x[15] = (kappa*phi2 - (phi2 - 1)*phi3 - 1)/norm - ratio*(-kappa*pi_c - pi_g - pi_a);
		}
		// dQdphi2
		else if(index == 2){
			const double dnorm = 2*kappa*phi1*phi1 + 4*(phi1*phi1 - (phi1*phi1 - 2*phi1 + 1)*phi2 - 2*phi1 + 1)*phi3*phi3 - 4*kappa*phi1 - 4*(kappa*phi1*phi1 - 2*kappa*phi1 + kappa)*phi2 + 2*((kappa + 3)*phi1 - 2*phi1*phi1 + 2*(kappa*phi1*phi1 - 2*kappa*phi1 + kappa)*phi2 - kappa - 1)*phi3 + 2*kappa;
			double ratio = dnorm/norm2;
			x[0] = (-kappa*(phi1 - 1)*phi3 + ((phi1 - 1)*phi3 - phi1 + 1) + (phi1 - 1))/norm - ratio*(-pi_t - kappa*pi_g - pi_c);
			x[1] = -(phi1 - 1)/norm - ratio*pi_c;
			x[2] = kappa*(phi1 - 1)*phi3/norm - ratio*kappa*pi_g;
			x[3] = -((phi1 - 1)*phi3 - phi1 + 1)/norm - ratio*pi_t;
			
			x[4] = -ratio*pi_a;
			x[5] = (-(phi1 - 1)*phi3 + ((phi1 - 1)*phi3 - phi1 + 1)*kappa)/norm - ratio*(-kappa*pi_t - pi_g - pi_a);
			x[6] = (phi1 - 1)*phi3/norm - ratio*pi_g;
			x[7] = -((phi1 - 1)*phi3 - phi1 + 1)*kappa/norm - ratio*kappa*pi_t;
			
			x[8] = -ratio*kappa*pi_a;
			x[9] = -(phi1 - 1)/norm - ratio*pi_c;
			x[10] = ((phi1 - 1) + (phi1 - 1)*phi3 - phi1 + 1)/norm - ratio*(-kappa*pi_a - pi_c - pi_t);
			x[11] = (-(phi1 - 1)*phi3 + phi1 - 1)/norm - ratio*pi_t;
			
			x[12] = -ratio*pi_a;
			x[13] = -kappa*(phi1 - 1)/norm - ratio*kappa*pi_c;
			x[14] = (phi1 - 1)*phi3/norm - ratio*pi_g;
			x[15] = (kappa*(phi1 - 1) - (phi1 - 1)*phi3)/norm - ratio*(-kappa*pi_c - pi_g - pi_a);
			
		}// dQdphi3
		else if(index == 3){
			const double dnorm = -2*(kappa - 2)*phi1*phi1 + 2*(kappa*phi1*phi1 - 2*kappa*phi1 + kappa)*phi2*phi2 + 2*(kappa - 3)*phi1 + 2*((kappa + 3)*phi1 - 2*phi1*phi1 - kappa - 1)*phi2 - 4*((phi1*phi1 - 2*phi1 + 1)*phi2*phi2 + phi1*phi1 - 2*(phi1*phi1 - 2*phi1 + 1)*phi2 - 2*phi1 + 1)*phi3 + 2;
			const double ratio = dnorm/norm2;
			x[0] = (-((phi1 - 1)*phi2 - phi1 + 1)*kappa + ((phi1 - 1)*phi2 - phi1 + 1))/norm - ratio*(-pi_t - kappa*pi_g - pi_c);
			x[1] = - ratio*pi_c;
			x[2] = ((phi1 - 1)*phi2 - phi1 + 1)*kappa/norm - ratio*kappa*pi_g;
			x[3] = -((phi1 - 1)*phi2 - phi1 + 1)/norm - ratio*pi_t;
			
			x[4] = - ratio*pi_a;
			x[5] = (-((phi1 - 1)*phi2 - phi1 + 1) + ((phi1 - 1)*phi2 - phi1 + 1)*kappa)/norm - ratio*(-kappa*pi_t - pi_g - pi_a);
			x[6] = ((phi1 - 1)*phi2 - phi1 + 1)/norm - ratio*pi_g;
			x[7] = -((phi1 - 1)*phi2 - phi1 + 1)*kappa/norm - ratio*kappa*pi_t;
			
			x[8] = - ratio*kappa*pi_a;
			x[9] = - ratio*pi_c;
			x[10] = ((phi1 - 1)*phi2 - phi1 + 1)/norm - ratio*(-kappa*pi_a - pi_c - pi_t);
			x[11] = (-(phi1 - 1)*phi2 + phi1 - 1)/norm - ratio*pi_t;
			
			x[12] = - ratio*pi_a;
			x[13] = - ratio*kappa*pi_c;
			x[14] = ((phi1 - 1)*phi2 - phi1 + 1)/norm - ratio*pi_g;
			x[15] = (-(phi1 - 1)*phi2 + phi1 - 1)/norm - ratio*(-kappa*pi_c - pi_g - pi_a);
		}
		else{
			exit(1);
		}
		m->dQ_need_update = false;
	}
	
	double *v = m->eigendcmp->eval;
	double xx[16];
	Matrix_mult3(mat, (const double**)m->eigendcmp->Invevec, x, 4,4,4,4);
	Matrix_mult4(xx, mat, (const double**)m->eigendcmp->evec, 4,4,4,4);
	// up to now the above operations can be recycled across branches
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			if(v[i] != v[j]){
				mat[i*4+j] = xx[i*4+j]*(exp(v[i]*t) - exp(v[j]*t))/(v[i]-v[j]);
			}
			else{
				mat[i*4+j] = xx[i*4+j]*t*exp(v[i]*t);
			}
		}
	}
	Matrix_mult3(xx, (const double**)m->eigendcmp->evec, mat, 4,4,4,4);
	Matrix_mult4(mat, xx, (const double**)m->eigendcmp->Invevec, 4,4,4,4);
}
