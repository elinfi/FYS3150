#ifndef GAUSS_QUADRATURE_H
#define GAUSS_QUADRATURE_H

double function(double x1, double y1, double z1, double x2, double y2, double z2);
double gauss_quadrature(int n, double lambda, bool timing);

#endif // GAUSS_QUADRATURE_H
