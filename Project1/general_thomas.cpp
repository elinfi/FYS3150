#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "time.h"
#include "functions.h"

using namespace std;

void gaussian_elimination (double* a, double* b, double* b_, double* c,
                             double* g, double* g_, double* v, int N) {

   // forward substitution
   b_[1] = b[1];
   g_[1] = g[1];
   for (int i=2; i<N+1; i++) {
       b_[i] = b[i] - (a[i-1]*c[i-1])/b_[i-1];
       g_[i] = g[i] - (g_[i-1]*a[i-1])/b_[i-1];
   }

   // backward substitution
   v[N] = g_[N]/b_[N];
   for (int i=N-1; i>0; i--) {
       v[i] = (g_[i] - a[i]*v[i+1])/b_[i];
   }
}

void general_thomas (int N) {

    // diagonals
    double* a = new double[N+2] ();      // below diagonal
    double* b = new double[N+2] ();      // leading diagonal
    double* b_ = new double[N+2] ();     // changed leading diagonal
    double* c = new double[N+2] ();      // above diagonal

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
        x[i] = i*h;
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;

        g[i] = hh*f(x[i]);
        u[i] = analytic_solution(x[i]);
    }

    // start timing
    auto start = chrono::high_resolution_clock::now();

    gaussian_elimination(a, b, b_, c, g, g_, v, N);

    // end timing
    auto finish = chrono::high_resolution_clock::now();

    // print time
    chrono::duration<double> elapsed = (finish - start);
    cout << "General Thomas Alg: " << elapsed.count() << " s\n";


    //Calculate relative error
    double* eps = new double[N];

    for (int i=1; i<N; i++){
        eps[i-1] = log10(abs((v[i]-u[i])/u[i]));
    }

    //Find maximum element in error array
    double max = eps[0];
    for (int i=1; i<N; i++){
        if (abs(eps[i]) > abs(max)) {
            max = eps[i];
        }
    }

    //Print max relative error
    cout << "N= "<< N << endl;
    cout << "log10(h): "<< log10(h) <<endl;
    cout << "Max rel. err.: " << max <<endl;


    // write the numerical solution to file
    string filename = "1b-";
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

    // write the analytic solution to file
    string file = "analytic_solution-";
    file.append(to_string(N));
    file.append(".txt");

    ofile.open(file);

    for (int i=0; i<N+2; i++) {
        if (i == N+1) {
            ofile << u[i];
        }
        else {
        ofile << u[i] << ",";
        }
    }
    ofile.close();

    delete[] a;
    delete[] b;
    delete[] b_;
    delete[] c;
    delete[] g;
    delete[] g_;
    delete[] v;
    delete[] u;
    delete[] x;
    delete[] eps;
}
