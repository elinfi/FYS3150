#include "monte_carlo.h"
#include <iostream>
#include <random>
#include <cmath>
#include <chrono>
#include <tuple>

using namespace std;

double int_function(double x1, double y1, double z1, double x2, double y2, double z2) {
    double value;
    double eps = 1.0E-10;

    double r1_r2 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
    if (r1_r2 < eps) {
        value = 0;
    }
    else {
        value = exp(-4*(sqrt(x1*x1 + y1*y1 + z1*z1) + sqrt(x2*x2 + y2*y2 + z2*z2)))*(1/r1_r2);
    }
    return value;
}


tuple<double, double> monte_carlo(int n, double lambda, int seed) {
    double a = -lambda;
    double b = lambda;

    // generate engine
    mt19937_64 engine(seed);

    // generate uniform distribution
    uniform_real_distribution<double> uniform(a, b);

    double sum = 0;
    double sum_sigma = 0;

    for (int i=0; i<n; i++) {
        double x1 = uniform(engine);
        double y1 = uniform(engine);
        double z1 = uniform(engine);
        double x2 = uniform(engine);
        double y2 = uniform(engine);
        double z2 = uniform(engine);

        double f = int_function(x1, y1, z1, x2, y2, z2);
        sum += f;
        sum_sigma += f*f;
    }

    double jacobi = pow(b-a, 6);
    sum = sum/n;

    sum_sigma  = sum_sigma/n;

    double variance = jacobi*(sum_sigma - sum*sum);

    return make_tuple(jacobi*sum, variance);
}
