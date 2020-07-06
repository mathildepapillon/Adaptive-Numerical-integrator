#ifndef NUM_INTEG_H
#define NUM_INTEG_H

double NumIntegrate(double *f, double h, int down, int up, int nmax); 

double NumIntegrateV(double (*func)(double), double x_down, double x_up, int num_intervals);

double AdaptiveIntegrate(double (*func)(double), double down, double h, double *f_prev, double ans_prev, double tol, int *count);

double Integrate(double (*integrand)(double), double x_down, double x_up, double tol, int *count, int limit);

#endif