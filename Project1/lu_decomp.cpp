#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "time.h"
#include "functions.h"

using namespace std;
using namespace arma;

void lu_decomposition (int N) {

    Mat<double> A = zeros(N,N);

    vec x(N);
    vec g(N);

    double h = 1.0/(N + 1);
    double hh = h*h;


    A(0, 0) = 2;
    for (int i = 0; i < N; i++) {
            if (i > 0) {
                A(i, i-1) = -1;        // immediately below leading diagonal
                A(i, i)   = 2;         // leading diagonal
            }
            if (i < N-1) {
                A(i, i+1) = -1;        // immediately above leading diagonal
            }
    }

    for (int i=1; i<N; ++i) {
        x[i] = i*h;
        g[i] = hh*f(x[i]);
    }

    mat L, U;

    // start timing
    auto start = chrono::high_resolution_clock::now();

    lu(L, U, A);

    vec y = solve(L, g);
    vec v = solve(U, y);

    // end timing
    auto finish = chrono::high_resolution_clock::now();

    // print time
    chrono::duration<double> elapsed = (finish - start);
    cout << "LU-decomposition: " << elapsed.count() << " s\n";

}
