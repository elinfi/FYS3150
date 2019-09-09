#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <cmath>
#include "time.h"
#include "functions.h"

using namespace std;


void optimized_gaussian_elimination (double* b_, double* g, double* g_, double* v, int N)
    {
   // forward substitution
   g_[1] = g[1];
   for (int i=2; i<N+1; i++) {
       g_[i] = g[i] + g_[i-1]/b_[i-1];
   }

   // backward substitution
   v[N] = g_[N]/b_[N];
   for (int i=N-1; i>0; i--) {
       v[i] = (g_[i] + v[i+1])/b_[i];
   }
}

void optimized_thomas (int N) {

    // diagonals
    double* b_ = new double[N+2] ();     // changed leading diagonal

    // right hand side
    double* g = new double[N+2] ();      // right side of equation
    double* g_ = new double[N+2] ();     // changed right side of equation

    // numerical solution
    double* v = new double[N+2] ();      // numerical solution

    // analytic solution
    double* u = new double[N+2] ();      // analytic solution


    // step size
    double h = 1.0/(N + 1);
    double hh = h*h;

    // steps
    double* x = new double[N+2] ();

    // initializing vectors
    for (int i=0; i < N+2; i++) {
        b_[i] = (i + 1.0)/i;

        x[i] = i*h;

        g[i] = hh*f(x[i]);
        u[i] = analytic_solution(x[i]);
    }


    // start timing
    auto start = chrono::high_resolution_clock::now();

    optimized_gaussian_elimination(b_, g, g_, v, N);

    // end timing
    auto finish = chrono::high_resolution_clock::now();

    // print time
    chrono::duration<double> elapsed = (finish - start);
    cout << "Optimized Thomas Alg: " << elapsed.count() << " s\n";


    // write the numerical solution to file
    string filename = "1c-";
    filename.append(to_string(N));
    filename.append(".txt");

    ofstream ofile;
    ofile.open(filename);

    for (int i=0; i<N+2; i++) {
        if (i == N+1) {
            ofile << v[i];
        }
        else {
        ofile << v[i] << ",";
        }
    }

    ofile.close();


    delete[] b_;
    delete[] g;
    delete[] g_;
    delete[] v;
    delete[] u;
    delete[] x;
}
