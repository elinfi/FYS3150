#ifndef IMPROVED_MONTE_CARLO_H
#define IMPROVED_MONTE_CARLO_H
#include <tuple>

using namespace std;


double improved_function(double r1, double r2, double theta1, double theta2, double phi1, double phi2);
tuple<double, double> improved_monte_carlo(int n, int seed);

#endif // IMPROVED_MONTE_CARLO_H
