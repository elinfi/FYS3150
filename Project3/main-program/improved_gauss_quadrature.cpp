#include "improved_gauss_quadrature.h"
#include "gauss_laguerre.h"
#include "gauss_legendre.h"
#include <cmath>
#include <iostream>
#include <chrono>

using namespace std;

double function_improved(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
    double cos_beta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    double root = r1*r1 + r2*r2 - 2*r1*r2*cos_beta;

    double eps = 10e-10;
    if (root < eps) {
        return 0;
    }

    double r12 = sqrt(root);

    double value = sin(theta1)*sin(theta2)/(1024*r12);

    return value;
}


double improved_gauss_quadrature(int n, bool timing) {
    double *r = new double [n+1];
    double *w_r = new double [n+1];

    double *theta = new double [n];
    double *phi = new double [n];
    double *w_theta = new double [n];
    double *w_phi = new double [n];

    double pi = acos(-1);

    // calculate the polynomials and weights
    gauss_laguerre(r, w_r, n, 2);
    gauss_legendre(0, pi, theta, w_theta, n);
    gauss_legendre(0, 2*pi, phi, w_phi, n);

    // start timing
    auto start = chrono::high_resolution_clock::now();

    // six for loops, one for respectively r1, r2, theta1, theta2, phi1, phi2
    double sum = 0;
    for (int i=1; i<n+1; i++) {
        for (int j=1; j<n+1; j++) {
            for (int k=0; k<n; k++) {
                for (int l=0; l<n; l++) {
                    for (int m=0; m<n; m++) {
                        for (int o=0; o<n; o++) {
                            sum += w_r[i]*w_r[j]*w_theta[k]*w_theta[l]*w_phi[m]*w_phi[o]*function_improved(r[i], r[j], theta[k], theta[l], phi[m], phi[o]);
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
        cout << "Improved Gauss Quadrature: " << elapsed.count() << " s\n";
    }


    delete [] r;
    delete [] w_r;
    delete [] theta;
    delete [] phi;
    delete [] w_theta;
    delete [] w_phi;

    return sum;


}
