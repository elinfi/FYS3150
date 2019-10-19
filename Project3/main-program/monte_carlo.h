#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include <tuple>

using namespace std;

double function(double x1, double y1, double z1, double x2, double y2, double z2);
tuple<double, double> monte_carlo(int n, double lambda, int seed);

#endif // MONTE_CARLO_H
