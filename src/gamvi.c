//
//  gamvi.c
//  viphy
//
//  Created by Mathieu Fourment on 12/04/2017.
//  Copyright © 2017 University of Technology Sydney. All rights reserved.
//

#include "gamvi.h"

#include <math.h>

#include "phyc/matrix.h"
#include "phyc/mstring.h"
#include "phyc/random.h"
#include "phyc/solve.h"


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>



//             for (int k = 0; k < dim; k++) {
//                 params[k+dim] = 5;
//                 params[k] = exp(Node_distance(nodes[k]))/params[k+dim];
//             }
//             
//             //for (int i = 0; i < paramCount; i++) params[i] = 0.1;
//             stochastic_gradient_ascent(tlk, nodes, dim, 1, 1, params, paramCount, eta, &elbo_gamma_meanfield, &grad_elbo_gamma_meanfield, r);
//             
//             printf("\n");
//             for( int i = 0; i < dim; i++ ){
//                 printf("%s t: %e a: %e b: %e mean: %e mode: %e\n", Node_name(nodes[i]), bls[i], exp(params[i]), exp(params[dim+i]), exp(params[i])/exp(params[i+dim]), (exp(params[i])-1)/exp(params[i+dim]));
//             }
//             printf("\n");
        
void grad_elbo_gamma_meanfield(SingleTreeLikelihood* tlk, Node **nodes, int nodeCount, int sampleCount, const double* params, double* grads,
                        const gsl_rng* r){
    double lambda = 20;
    int dim = nodeCount;

    memset(grads, 0, sizeof(double)*dim*2);
    
    double *dlnlt = dvector(Tree_node_count(tlk->tree));
    double eps = 1.0e-10;
    
    double* z = dvector(dim);
    double* u = dvector(dim);
    bool reparam = true;
    
    for (int i = 0; i < sampleCount; i++) {
        for (int k = 0; k < dim; k++) {
            
            double alpha = exp(params[k]);
            double beta = exp(params[k+dim]);
            
            if(reparam){
                // gsl uses shape and scale to parametrize gamma distribution
                z[k] = gsl_ran_gamma(r, alpha, 1.0/beta);
            }
            else{
                // a,b > 0
                u[k] = gsl_rng_uniform(r);
                
                // psi == logarithmetic derivative of gamma function
                if(alpha < 1 && (24-22.6*u[k])*log(alpha) < -10){
                    z[k] = exp((log(u[k])+log(alpha) + gsl_sf_lngamma(alpha))/alpha - log(beta)); // F^{-1}_{a,b}(z)=x
                }
                else{
                    z[k] = gsl_cdf_gamma_Pinv(u[k], alpha, 1.0/beta);
                }
            }
            Node_set_distance(nodes[k], z[k]);
            SingleTreeLikelihood_update_one_node(tlk, nodes[k]);
        }
        
        // Tree log(likelihood * prior)
        double logP = tlk->calculate(tlk);
        logP += log(lambda)*dim;
        for (int k = 0; k < dim; k++) {
            logP -= lambda*z[k];
        }
        
        // Calculate 1st derivative for each branch
        calculate_all_dlnl_dt(tlk, dlnlt);
        
        for (int k = 0; k < dim; k++){
            Node* theNode = nodes[k];
            
            const double dprior = -lambda; // first derivative of log(prior) = prior'/prior
            double dlnl = dlnlt[Node_id(theNode)]; // first derivative of loglikelihood function
            double dlogP = dlnl + dprior;
            
            double alpha = exp(params[k]);
            double beta = exp(params[k+dim]);
            
            if(reparam){
                double psi_alpha = gsl_sf_psi(alpha);
                double psi1_alpha = gsl_sf_psi_1(alpha);
                double psi2_alpha = gsl_sf_psi_n(2, alpha);
                
                double epsilon = (log(z[k]) - psi_alpha + log(beta))/sqrt(psi1_alpha);
                
                double h_alpha = z[k] * (epsilon*psi2_alpha/(2.0*sqrt(psi1_alpha)) + psi1_alpha);
                double u_alpha = epsilon*psi2_alpha/(2.0*sqrt(psi1_alpha)) + psi1_alpha + (psi2_alpha/(2.0*psi1_alpha));
                
                double dlogqdz = (alpha-1.0)/z[k] - beta;
                
                // alpha
                double g_rep = dlogP * h_alpha;
                double dlogqda = log(beta) - psi_alpha + log(z[k]);
                double g_corr = logP*(dlogqdz*h_alpha + dlogqda + u_alpha);
                grads[k] += g_rep + g_corr;
                
                // beta
                // g_corr == 0
                // g_rep
                grads[k+dim] += -dlogP*z[k]/beta;
                
            }
            else{
                double dfda;
                double dfdb = -z[k]/beta;
                if(alpha < 1 && (24-22.6*u[k])*log(alpha) < -10){
                    double dlogxda = -(log(z[k]) + log(alpha) + gsl_sf_lngamma(alpha))/(alpha*alpha) + (1.0/alpha + gsl_sf_psi(alpha))/alpha;
                    dfda = dlogxda * z[k];
                }
                else{
                    dfda = (gsl_cdf_gamma_Pinv(u[k], alpha+eps, 1.0/beta) - z[k])/eps;
                }
                dlogP -= (alpha-1.0)/z[i] - beta;
                grads[k] += dlogP * dfda;
                grads[k+dim] += dlogP * dfdb;
            }
        }
    }
    
    for (int k = 0; k < dim*2; k++) {
        grads[k] /= sampleCount;
    }
    
    if(reparam){
        // Gradient entropy
        for (int k = 0; k < dim; k++) {
            double alpha = exp(params[k]);
            double beta = exp(params[k+dim]);
            grads[k] += 1.0 + (1.0-alpha) * gsl_sf_psi_1(alpha);
            grads[k+dim] += -1.0/beta;
        }
    }
    
    free(dlnlt);
    free(z);
    free(u);
}

double elbo_gamma_meanfield(SingleTreeLikelihood* tlk, Node **nodes, int nodeCount, int sampleCount, const double* params, const gsl_rng* r){
    double lambda = 20;
    int dim = nodeCount;
    double elbo = 0;
    bool reparam = true;
    
    for (int i = 0; i < sampleCount; i++) {
        for (int k = 0; k < dim; k++) {
            
            double alpha = exp(params[k]);
            double beta = exp(params[k+dim]);
            
            if(reparam){
                double z = gsl_ran_gamma(r, alpha, 1.0/beta);
                Node_set_distance(nodes[k], z );
                SingleTreeLikelihood_update_one_node(tlk, nodes[k]);
            }
            else{
                double u = gsl_rng_uniform(r);
                double z;
                
                // psi == logarithmetic derivative of gamma function
                if(alpha < 1 && (24-22.6*z)*log(alpha) < -10){
                    z = exp((log(u)+log(alpha) + gsl_sf_lngamma(alpha))/alpha - log(beta)); // F^{-1}_{a,b}(z)=x
                }
                else{
                    z = gsl_cdf_gamma_Pinv(u, alpha, 1.0/beta);
                }
                Node_set_distance(nodes[k], z);
                SingleTreeLikelihood_update_one_node(tlk, nodes[k]);
            }
        }
        
        elbo += tlk->calculate(tlk);
        // Calculate prior
        elbo += log(lambda)*dim;
        for (int k = 0; k < dim; k++) {
            elbo -= lambda*Node_distance(nodes[k]);
        }
    }
    elbo /= sampleCount;
    
    // Entropy of gamma: H(X) = alpha - ln(beta) lngamma(alpha) + (1-a)psi(alpha)
    for (int k = 0; k < dim; k++) {
        double alpha = exp(params[k]);
        elbo += alpha - params[k+dim] + gsl_sf_lngamma(alpha) + (1.0-alpha) * gsl_sf_psi(alpha);
    }
    return elbo;
}


void grad_elbo_gamma_meanfield2(SingleTreeLikelihood* tlk, Node **nodes, int nodeCount, int sampleCount, const double* params, double* grads,
                                const gsl_rng* r){
    double lambda = 20;
    int dim = nodeCount;
    memset(grads, 0, sizeof(double)*dim*2);
    
    
    for (int i = 0; i < sampleCount; i++) {
        
        double logQ = 0;
        
        for (int k = 0; k < dim; k++) {
            // k lambda alpha <-> T(alpha)
            // k+dim lambda mu <-> T(mu)/T(alpha)
            double alpha = softplus(params[k]);
            double mu = softplus(params[k+dim]);
            //double theta = gsl_ran_gamma(r, alpha, mu/alpha);
            double theta = gsl_ran_gamma(r, alpha, alpha/mu);
            Node_set_distance(nodes[k], theta);
            SingleTreeLikelihood_update_one_node(tlk, nodes[k]);
            logQ += alpha * log(alpha) - alpha * log(mu) - gsl_sf_lngamma(alpha) + (alpha-1.) * log(theta) - alpha * theta / mu;//log(gsl_ran_gamma_pdf(theta, alpha, mu/alpha));
        }
        //Tree_print_newick(stdout, tlk->tree, false);
        double logP = tlk->calculate(tlk);
        // Calculate prior
        double prior = log(lambda)*dim;
        for (int k = 0; k < dim; k++) {
            Node* node = nodes[k];
            if(!Node_isroot(node) && Node_right(Tree_root(tlk->tree)) != node){
                prior -= lambda*Node_distance(node);
            }
        }
        
        logP += prior;
        
        double p_q = logP - logQ;
        
        for (int k = 0; k < dim; k++){
            Node* theNode = nodes[k];
            
            double x = Node_distance(theNode);
            
            double alpha = softplus(params[k]);
            double mu = softplus(params[k+dim]);
            double grad_alpha = log(alpha) + 1.0 - log(mu) - gsl_sf_psi(alpha) + log(x) - (x/mu);
            double grad_mu = -(alpha/mu)-(alpha*x)/(mu*mu);
            double grad_lambda_alpha = exp(params[k])/(1.0+exp(params[k])) * grad_alpha;
            double grad_lambda_mu = exp(params[k+dim])/(1.0+exp(params[k+dim])) * grad_mu;
            
            grads[k] += grad_lambda_alpha * p_q;
            grads[dim+k] += grad_lambda_mu * p_q;
        }
    }
    
    for (int k = 0; k < dim*2; k++) {
        grads[k] /= sampleCount;
    }
    
    //print_dvector(grads, dim*2);
    
}

double elbo_gamma_meanfield2(SingleTreeLikelihood* tlk, Node **nodes, int nodeCount, int sampleCount, const double* params,
                            const gsl_rng* r){
    double lambda = 20;
    int dim = nodeCount;
    double elbo = 0;
    
    for (int i = 0; i < sampleCount; i++) {
        
        double logQ = 0;
        
        for (int k = 0; k < dim; k++) {
            // k lambda alpha <-> T(alpha)
            // k+dim lambda mu <-> T(mu)/T(alpha)
            double alpha = softplus(params[k]);
            double mu = softplus(params[k+dim]);
            //double theta = gsl_ran_gamma(r, alpha, mu/alpha);
            double theta = gsl_ran_gamma(r, alpha, alpha/mu);
            Node_set_distance(nodes[k], theta);
            SingleTreeLikelihood_update_one_node(tlk, nodes[k]);
            logQ += alpha * log(alpha) - alpha * log(mu) - gsl_sf_lngamma(alpha) + (alpha-1.) * log(theta) - alpha * theta / mu;//log(gsl_ran_gamma_pdf(theta, alpha, mu/alpha));
        }
        Tree_print_newick(stdout, tlk->tree, false);
        double logP = tlk->calculate(tlk);
        // Calculate prior
        double prior = log(lambda)*dim;
        for (int k = 0; k < dim; k++) {
            Node* node = nodes[k];
            if(!Node_isroot(node) && Node_right(Tree_root(tlk->tree)) != node){
                prior -= lambda*Node_distance(node);
            }
        }
        
        logP += prior;
        
        elbo += (logP - logQ)*exp(logQ);
    }
    return elbo/sampleCount;
}

void stochastic_gradient_ascentRMSprop(SingleTreeLikelihood* tlk, Node **nodes, int nodeCount, int elbo_samples, int grad_samples,
                                       double *parameters, int parameterCount,
                                       double eta,
                                       double(*elbofn)(SingleTreeLikelihood*, Node**, int, int, const double*),
                                       void(*grad_elbofn)(SingleTreeLikelihood*, Node**, int, int, const double*, double*),
                                       const gsl_rng* r){
    
    double *grads = calloc(parameterCount, sizeof(double));
    double *history_grad_squared = calloc(parameterCount, sizeof(double));
    
    int iter = 1;
    int max_iter = 10000;
    double tau = 1;
    double pre_factor  = 0.9;
    double post_factor = 0.1;
    double elbo = 0;
    double elbo_prev = -INFINITY;
    double elbo_best = -INFINITY;
    double tol_rel_obj = 0.001;
    int eval_elbo = 10;
    
    double *elbos = calloc(max_iter/eval_elbo, sizeof(double));
    
    while(iter++ < max_iter){
        grad_elbofn(tlk, nodes, nodeCount, grad_samples, parameters, grads);
        
        double eta_scaled = eta / sqrt(iter);
        
        // Update step-size
        if (iter == 1) {
            for (int i = 0; i < parameterCount; i++) {
                history_grad_squared[i] = grads[i]*grads[i];
            }
        } else {
            for (int i = 0; i < parameterCount; i++) {
                history_grad_squared[i] = pre_factor * history_grad_squared[i] + post_factor * grads[i]*grads[i];
            }
        }
        double rho = pow((iter + 1024), -0.7);
        for (int i = 0; i < parameterCount; i++) {
            parameters[i] += rho * grads[i] / (sqrt(history_grad_squared[i]));
        }
        print_dvector(parameters, parameterCount);
        
        //if (iter % eval_elbo == 0) {
        elbo_prev = elbo;
        elbo = elbofn(tlk, nodes, nodeCount, elbo_samples, parameters);
        printf("%d ELBO: %f (%f)\n",iter, elbo, elbo_prev);
        
        //        if (elbo > elbo_best){
        //            elbo_best = elbo;
        //        }
        //        double delta_elbo = fabs((elbo_prev - elbo) / elbo);
        //        elbos[iter/eval_elbo-1] = delta_elbo;
        //        qsort (elbos, iter/eval_elbo, sizeof(double), compare);
        //        int m = iter/eval_elbo/2;
        //        if(elbos[m] < tol_rel_obj){
        //            printf("ELBO converged: %f < %f  %f  %d %d\n",elbos[m], tol_rel_obj, delta_elbo, m, iter);
        //            //break;
        //        }
        //}
    }
    free(elbos);
    free(history_grad_squared);
    free(grads);
}
