#ifndef HEADER_H
#define HEADER_H

#endif // HEADER_H

double legendre_poly(int n, double x);
void gauleg(double x1, double x2, double x[], double w[], int n);
double int_func(double x1, double x2, double y1, double y2, double z1, double z2);
void gauss_laguerre(double *x, double *w, int n, double alf);
double laguerre_func(double r1, double r2, double thet1, double thet2, double phi1, double phi2);
double montecarlo_func(double lambda, int n, double & variance);

void outfile_c(double res, double lambda, int n);
double improved_montecarlo_func(double r1_, double r2_, double thet1, double thet2, double phi1, double phi2);
double montecarlo_two_electric_boogaloo(int n, int seed, double & variance);

