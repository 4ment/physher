//
//  test_distributions.c
//  physher
//  Created by Mathieu Fourment on 13/03/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include <gsl/gsl_randist.h>

#include "minunit.h"

#include "phyc/distexp.h"
#include "phyc/distgamma.h"
#include "phyc/distlognormal.h"
#include "phyc/distnormal.h"
#include "phyc/parameters.h"

#pragma region Exponential Distribution

char* test_exponential_distribution_aux(Parameter* x, Parameter* lambda, double* dlogPdx, double* dlogPdlambda) {
    const double *lambdaValues = Parameter_values(lambda);
    const double *xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(1);
    Parameters_move(parameters, lambda);

    DistributionModel* dm = new_ExponentialDistributionModel_with_parameters(parameters, xs, DISTRIBUTION_EXPONENTIAL_RATE);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if(Parameter_size(lambda) > 1){
        for(size_t i = 0; i < Parameter_size(x); i++) {
            logP2 += log(gsl_ran_exponential_pdf(xValues[i], 1/lambdaValues[i]));
        }
    }
    else{
        for(size_t i = 0; i < 3; i++) {
            logP2 += log(gsl_ran_exponential_pdf(xValues[i], 1/lambdaValues[0]));
        }
    }
    mu_assert(logP2 == logP, "logP not matching");
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    Parameters*ps = new_Parameters(2);
    Parameters_add(ps, lambda);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for(size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for(size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for(size_t i = 0; i < Parameter_size(lambda); i++) {
        mu_assert(fabs(lambda->grad[i] - dlogPdlambda[i]) < 0.0001, "dlogPdlambda not matching");
    }

    model->free(model);
    return NULL;
}


char* test_exponential_distribution_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double lambdaValue = 2.0;
    Parameter* lambda = new_Parameter("lambda", lambdaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-2.0, -2.0, -2.0};
    double dlogPdlambda = -1.0099999999999998;

    return test_exponential_distribution_aux(x, lambda, dlogPdx, &dlogPdlambda);
}

char* test_exponential_distribution() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double lambdaValues[] = {0.1, 1.0, 2.0};
    Parameter* lambda = new_Parameter2("lambda", lambdaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-0.1, -1.0, -2.0};
    double dlogPdlambda[3] = {8.0, 0.5, 0.49};

    return test_exponential_distribution_aux(x, lambda, dlogPdx, dlogPdlambda);
}

#pragma endregion

#pragma region Normal Distribution

char* test_normal_distribution_aux(Parameter* x, Parameter* mu, Parameter* sigma, double* dlogPdx, double* dlogPdmu, double* dlogPdsigma) {
    const double *muValues = Parameter_values(mu);
    const double *sigmaValues = Parameter_values(sigma);
    const double *xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, mu);
    Parameters_move(parameters, sigma);

    DistributionModel* dm = new_NormalDistributionModel_with_parameters(parameters, xs, DISTRIBUTION_NORMAL_MEAN_SIGMA);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if(Parameter_size(mu) > 1){
        for(size_t i = 0; i < Parameter_size(x); i++) {
            logP2 += log(gsl_ran_gaussian_pdf(xValues[i] - muValues[i], sigmaValues[i]));
        }
    }
    else{
        for(size_t i = 0; i < 3; i++) {
            logP2 += log(gsl_ran_gaussian_pdf(xValues[i] - muValues[0], sigmaValues[0]));
        }
    }
    mu_assert(logP2 == logP, "logP not matching");
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    Parameters*ps = new_Parameters(3);
    Parameters_add(ps, mu);
    Parameters_add(ps, sigma);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for(size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for(size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for(size_t i = 0; i < Parameter_size(mu); i++) {
        mu_assert(fabs(mu->grad[i] - dlogPdmu[i]) < 0.0001, "dlogPdmu not matching");
        mu_assert(fabs(sigma->grad[i] - dlogPdsigma[i]) < 0.0001, "dlogPdsigma not matching");
    }

    model->free(model);
    return NULL;
}

char* test_normal_distribution_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double muValue = 2.0;
    double sigmaValue = 0.5;
    Parameter* mu = new_Parameter("mu", muValue, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma = new_Parameter("sigma", sigmaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-0.0, 6.0, 7.96};
    double dlogPdmu = -13.96;
    double dlogPdsigma = 43.680800000000005;

    return test_normal_distribution_aux(x, mu, sigma, dlogPdx, &dlogPdmu, &dlogPdsigma);
}

char* test_normal_distribution() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double muValues[] = {1.0, 2.0, 3.0};
    double sigmaValues[] = {0.5, 0.1, 0.2};
    Parameter* mu = new_Parameter2("mu", muValues, 3, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma = new_Parameter2("sigma", sigmaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-4.0, 149.99999999999997, 74.75};
    double dlogPdmu[3] = {4.0, -149.99999999999997, -74.75};
    double dlogPdsigma[3] = {6.0, 2239.999999999999, 1112.5124999999996};

    return test_normal_distribution_aux(x, mu, sigma, dlogPdx, dlogPdmu, dlogPdsigma);
}
#pragma endregion

#pragma region Lognormal Distribution

char* test_lognormal_distribution_aux(Parameter* x, Parameter* mu, Parameter* sigma, double* dlogPdx, double* dlogPdmu, double* dlogPdsigma) {
    const double *muValues = Parameter_values(mu);
    const double *sigmaValues = Parameter_values(sigma);
    const double *xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, mu);
    Parameters_move(parameters, sigma);

    DistributionModel* dm = new_LogNormalDistributionModel_with_parameters(parameters, xs, DISTRIBUTION_LOGNORMAL_MU_SIGMA);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if(Parameter_size(mu) > 1){
        for(size_t i = 0; i < Parameter_size(x); i++) {
            logP2 += log(gsl_ran_lognormal_pdf(xValues[i], muValues[i], sigmaValues[i]));
        }
    }
    else{
        for(size_t i = 0; i < 3; i++) {
            logP2 += log(gsl_ran_lognormal_pdf(xValues[i], muValues[0], sigmaValues[0]));
        }
    }
    mu_assert(logP2 == logP, "logP not matching");
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    Parameters*ps = new_Parameters(3);
    Parameters_add(ps, mu);
    Parameters_add(ps, sigma);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for(size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for(size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for(size_t i = 0; i < Parameter_size(mu); i++) {
        mu_assert(fabs(mu->grad[i] - dlogPdmu[i]) < 0.0001, "dlogPdmu not matching");
        mu_assert(fabs(sigma->grad[i] - dlogPdsigma[i]) < 0.0001, "dlogPdsigma not matching");
    }

    model->free(model);
    return NULL;
}

char* test_lognormal_distribution_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double muValue = 2.0;
    double sigmaValue = 0.5;
    Parameter* mu = new_Parameter("mu", muValue, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma = new_Parameter("sigma", sigmaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {2.113705638880109, 19.545177444479563, 2542.068074395236};
    double dlogPdmu = -42.42068074395236;
    double dlogPdsigma = 414.7134337096188;

    return test_lognormal_distribution_aux(x, mu, sigma, dlogPdx, &dlogPdmu, &dlogPdsigma);
}

char* test_lognormal_distribution() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double muValues[] = {1.0, 2.0, 3.0};
    double sigmaValues[] = {0.5, 0.1, 0.2};
    Parameter* mu = new_Parameter2("mu", muValues, 3, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma = new_Parameter2("sigma", sigmaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {0.11370563888010943, 536.629436111989, 18912.925464970223};
    double dlogPdmu[3] = {-1.2274112777602189, -269.3147180559945, -190.12925464970223};
    double dlogPdsigma[3] = {-1.2467307776135135, 7243.041736157981, 7224.826694730264};

    return test_lognormal_distribution_aux(x, mu, sigma, dlogPdx, dlogPdmu, dlogPdsigma);
}

char* test_lognormal_distribution_prior_mean_stdev() {
    Parameter* m = new_Parameter("mean", 2.0, new_Constraint(0, INFINITY));
    Parameter* stdev = new_Parameter("stdev", 0.5, new_Constraint(0, INFINITY));
    double mValue = Parameter_value(m);
    double stdevValue = Parameter_value(stdev);

    Parameter* x1 = new_Parameter("x1", 0.5, new_Constraint(0, INFINITY));
    double x3Values[] = {2.0, 0.5, 0.01};
    Parameter* x3 = new_Parameter2("x3", x3Values, 3, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x3);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, m);
    Parameters_move(parameters, stdev);

    DistributionModel* dm = new_LogNormalDistributionModel_with_parameters(parameters, xs, DISTRIBUTION_LOGNORMAL_MEAN_STDEV);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    for(size_t i = 0; i < 3; i++) {
        double sigma = sqrt(log(1.0 + (stdevValue * stdevValue) / (mValue * mValue)));
        double mu = log(mValue) - 0.5 * sigma * sigma;
        logP2 += log(gsl_ran_lognormal_pdf(x3Values[i], mu, sigma));
    }
    mu_assert(logP2 == logP, "logP not matching");
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    Parameters*ps = new_Parameters(3);
    Parameters_add(ps, m);
    Parameters_add(ps, stdev);
    Parameters_add(ps, x3);
    model->gradient(model, ps);
    double dlogPdx[3] = {-0.7499999999999998, 42.73370751301173, 8589.547081367036};
    double dlogPdm = -292.9289997760196;
    double dlogPdstdev = 954.1913499637261;

    mu_assert(fabs(x3->grad[0] - dlogPdx[0]) < 0.0001, "dlogPdx1 not matching");
    mu_assert(fabs(x3->grad[1] - dlogPdx[1]) < 0.0001, "dlogPdx2 not matching");
    mu_assert(fabs(x3->grad[2] - dlogPdx[2]) < 0.0001, "dlogPdx3 not matching");
    mu_assert(fabs(m->grad[0] - dlogPdm) < 0.0001, "dlogPdm not matching");
    mu_assert(fabs(stdev->grad[0] - dlogPdstdev) < 0.0001, "dlogPdstdev not matching");

    printf("Gradients:\n");
    for (size_t i = 0; i < 3; i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for(size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    model->free(model);
    return NULL;
}

#pragma endregion

#pragma region Gamma Distribution

char* test_gamma_distribution_aux(Parameter* x, Parameter* alpha, Parameter* beta, double* dlogPdx, double* dlogPdalpha, double* dlogPdbeta, distribution_parameterization parameterization) {
    const double *alphaValues = Parameter_values(alpha);
    const double *betaValues = Parameter_values(beta);
    const double *xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, alpha);
    Parameters_move(parameters, beta);

    DistributionModel* dm = new_GammaDistributionModel_with_parameters(parameters, xs, parameterization);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if(Parameter_size(alpha) > 1){
        for(size_t i = 0; i < Parameter_size(x); i++) {
            double betaValue = parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE ? 1.0/betaValues[i] : betaValues[i];
            logP2 += log(gsl_ran_gamma_pdf(xValues[i], alphaValues[i], betaValue));
        }
    }
    else{
        double betaValue = parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE ? 1.0/betaValues[0] : betaValues[0];
        for(size_t i = 0; i < 3; i++) {
            logP2 += log(gsl_ran_gamma_pdf(xValues[i], alphaValues[0], betaValue));
        }
    }
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    mu_assert(logP2 == logP, "logP not matching");
    Parameters*ps = new_Parameters(3);
    Parameters_add(ps, alpha);
    Parameters_add(ps, beta);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for(size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for(size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for(size_t i = 0; i < Parameter_size(alpha); i++) {
        mu_assert(fabs(alpha->grad[i] - dlogPdalpha[i]) < 0.0001, "dlogPdalpha not matching");
        mu_assert(fabs(beta->grad[i] - dlogPdbeta[i]) < 0.0001, "dlogPdbeta not matching");
    }

    model->free(model);
    return NULL;
}

char* test_gamma_distribution_shape_rate_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double alphaValue = 2.0;
    double betaValue = 0.5;
    Parameter* alpha = new_Parameter("alpha", alphaValue, new_Constraint(0, INFINITY));
    Parameter* beta = new_Parameter("beta", betaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {0.0, 1.5, 99.5};
    double dlogPdalpha = -7.952964732963327;
    double dlogPdbeta = 9.49;

    return test_gamma_distribution_aux(x, alpha, beta, dlogPdx, &dlogPdalpha, &dlogPdbeta, DISTRIBUTION_GAMMA_SHAPE_RATE);
}

char* test_gamma_distribution_shape_rate() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double alphaValues[] = {1.0, 2.0, 3.0};
    double betaValues[] = {0.5, 0.1, 0.2};
    Parameter* alpha = new_Parameter2("alpha", alphaValues, 3, new_Constraint(0, INFINITY));
    Parameter* beta = new_Parameter2("beta", betaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-0.5, 1.9, 199.8};
    double dlogPdalpha[3] = {0.5772156649015329, -3.4185166086524577, -7.137392433520659};
    double dlogPdbeta[3] = {0.0, 19.5, 14.99};
    return test_gamma_distribution_aux(x, alpha, beta, dlogPdx, dlogPdalpha, dlogPdbeta, DISTRIBUTION_GAMMA_SHAPE_RATE);
}

char* test_gamma_distribution_shape_scale() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double alphaValues[] = {1.0, 2.0, 3.0};
    double betaValues[] = {0.5, 0.1, 0.2};
    Parameter* alpha = new_Parameter2("alpha", alphaValues, 3, new_Constraint(0, INFINITY));
    Parameter* beta = new_Parameter2("beta", betaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-2.0, -8.0, 195.0};
    double dlogPdalpha[3] = {1.9635100260214235, 1.1866535773356337, -3.9185166086524577};
    double dlogPdbeta[3] = {6.0, 30.0, -14.75};
    return test_gamma_distribution_aux(x, alpha, beta, dlogPdx, dlogPdalpha, dlogPdbeta, DISTRIBUTION_GAMMA_SHAPE_SCALE);
}

char* test_gamma_distribution_shape_scale_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double alphaValue = 2.0;
    double betaValue = 0.5;
    Parameter* alpha = new_Parameter("alpha", alphaValue, new_Constraint(0, INFINITY));
    Parameter* beta = new_Parameter("beta", betaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-1.5, 0.0, 98.0};
    double dlogPdalpha = -3.794081649603656;
    double dlogPdbeta = -1.9600000000000009;

    return test_gamma_distribution_aux(x, alpha, beta, dlogPdx, &dlogPdalpha, &dlogPdbeta, DISTRIBUTION_GAMMA_SHAPE_SCALE);
}

#pragma endregion

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_exponential_distribution);
    mu_run_test(test_exponential_distribution_prior);
    
    mu_run_test(test_normal_distribution);
    mu_run_test(test_normal_distribution_prior);

    mu_run_test(test_lognormal_distribution);
    mu_run_test(test_lognormal_distribution_prior);
    mu_run_test(test_lognormal_distribution_prior_mean_stdev);

    mu_run_test(test_gamma_distribution_shape_rate);
    mu_run_test(test_gamma_distribution_shape_rate_prior);
    mu_run_test(test_gamma_distribution_shape_scale);
    mu_run_test(test_gamma_distribution_shape_scale_prior);

    return NULL;
}

RUN_TESTS(all_tests);
