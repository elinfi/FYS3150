#include "transactions.h"
#include <iostream>
#include <random>
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;


void transactions(int N, double m0, double txn, ofstream& ofile) {
    int seed = time(0);

    // generate engine
    mt19937_64 engine(seed);

    // generate integers in interval [0, N-1] with uniform distribution
    uniform_int_distribution<int> int_uniform(0, N-1);

    uniform_int_distribution<int> int2_uniform(0, N-2);

    // generate double with uniform distribution
    uniform_real_distribution<double> double_uniform(0, 1);

    // Store the 500 agents in a vector with initial amount of money m0
    vec m = zeros<vec>(N);
    m.fill(m0);

    for (int t=0; t<txn; t++) {
        // choose a random pair (i, j)
        int i = int_uniform(engine);
        int j = int2_uniform(engine);
        if (j >= i) {
            j += 1;
        }

        // choose a random reassignment eps
        double eps = double_uniform(engine);

        // after transaction
        double mi_new = eps*(m[i] + m[j]);
        double mj_new = (1 - eps)*(m[i] + m[j]);

        m[i] = mi_new;
        m[j] = mj_new;
    }

    for (int i=0; i<N-1; i++) {
        ofile << m[i] << " ";
    }
    ofile << "\n";

}
