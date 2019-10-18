#include "improved_monte_carlo.h"
#include <cmath>
#include <iostream>
#include <random>
#include <chrono>

using namespace std;

double improved_function(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
    double cos_beta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    double root = r1*r1 + r2*r2 - 2*r1*r2*cos_beta;

    double eps = 10e-10;
    if (root < eps) {
        return 0;
    }

    double r12 = sqrt(root);

    double value = r1*r1*r2*r2*sin(theta1)*sin(theta2)/(1024*r12);

    return value;
}

double improved_monte_carlo(int n, bool timing) {
    double pi = acos(-1);
    double a = 0;
    double b_theta = pi;
    double b_phi = 2*pi;

    // generate engine
//    int seed = time_t(0);
    int seed = 1424;
    mt19937_64 engine(seed);

    // generate uniform distribution
    uniform_real_distribution<double> uniform(0, 1);

    // generate exponential distribution
    exponential_distribution<double> exponential(1);


    double sum = 0;
    double sum_sigma = 0;

    // start timing
    auto start = chrono::high_resolution_clock::now();

    for (int i=0; i<n; i++) {
        double r1 = exponential(engine);
        double r2 = exponential(engine);

        // use the mapping z = a + (b - a)*x, where x lies between 0 and 1
        double theta1 = a + (b_theta - a)*uniform(engine);
        double theta2 = a + (b_theta - a)*uniform(engine);
        double phi1 = a + (b_phi - a)*uniform(engine);
        double phi2 = a + (b_phi - a)*uniform(engine);

        double f = improved_function(r1, r2, theta1, theta2, phi1, phi2);

        sum += f;
        sum_sigma += f*f;
    }
    double jacobi = pow(b_theta - a, 2)*pow(b_phi - a, 2);
    sum = sum/n;

    // end timing
    auto finish = chrono::high_resolution_clock::now();

    if (timing) {
        // print time
        chrono::duration<double> elapsed = (finish - start);
        cout << "Improved Monte Carlo: " << elapsed.count() << " s\n";
    }


    sum_sigma  = sum_sigma/n;

    double variance = jacobi*(sum_sigma - sum*sum);

    cout << "variance: " << variance << endl;

    return jacobi*sum;
}
