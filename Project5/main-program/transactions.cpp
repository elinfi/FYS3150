#include "transactions.h"
#include <iostream>
#include <random>
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;


void transactions(int N, double m0, double lambda, double txn, ofstream& ofile) {
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
    vec m_mean = zeros<vec>(N);
    m_mean.fill(m0);

    for (int n=0; n<N; n++) {
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
            double dm = (1 - lambda)*(eps*m[j] - (1 - eps)*m[i]);
            m[i] += dm;
            m[j] -= dm;

            m_mean[i] += m[i];
            m_mean[j] += m[j];
        }
    }

    m_mean /= N;

    double mean = 0;
    double mean2 = 0;
    for (int i=0; i<N-1; i++) {
        mean += m_mean[i];
        mean2 += m_mean[i]*m_mean[i];
//        ofile << m_mean[i] << " ";
    }
//    ofile << "\n";

    mean /= N;
    mean2 /= N;

    double var = mean2 - mean*mean;
    ofile <<txn << " " << var << endl;

}
