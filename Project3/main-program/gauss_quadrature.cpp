#include "gauss_quadrature.h"
#include "gauss_legendre.h"
#include <iostream>
#include <cmath>
#include <chrono>

using namespace std;

double function(double x1, double y1, double z1, double x2, double y2, double z2) {
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

double gauss_quadrature(int n, double lambda, bool timing) {
    double *x = new double [n];
    double *w = new double [n];

    gauss_legendre(-lambda, lambda, x, w, n);

    // start timing
    auto start = chrono::high_resolution_clock::now();

    // six for loops, one for respectively x1, y1, z1, x2, y2, z2
    double sum = 0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            for (int k=0; k<n; k++) {
                for (int l=0; l<n; l++) {
                    for (int m=0; m<n; m++) {
                        for (int o=0; o<n; o++) {
                            sum += w[i]*w[j]*w[k]*w[l]*w[m]*w[o]*function(x[i], x[j], x[k], x[l], x[m], x[o]);
                        }
                    }
                }
            }
        }
    }

    // end timing
    auto finish = chrono::high_resolution_clock::now();

    if (timing) {
        // print time
        chrono::duration<double> elapsed = (finish - start);
        cout << "Gauss Quadrature: " << elapsed.count() << " s\n";
    }

    delete [] x;
    delete [] w;

    return sum;
}
