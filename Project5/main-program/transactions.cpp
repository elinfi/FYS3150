#include "transactions.h"
#include <iostream>
#include <random>
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

void test_transactions(int N, double m0, vec &m);

void transactions(int N, double m0, double lambda, double alpha, double gamma, double txn, ofstream& ofile) {
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

    // Matrix to store the number of interactions between two agents
    Mat<double> C(N, N, fill::zeros);

    for (int t=0; t<txn; t++) {
        // choose a random pair (i, j)
        int i = int_uniform(engine);
        int j = int2_uniform(engine);
        if (j >= i) {
            j += 1;
        }

        double p;
        if (m[i] == m[j]) {
            p = pow(C(min(i, j), max(i, j)) + 1, gamma);
        }
        else {
            p = pow(abs(m[i] - m[j]), -alpha)*pow(C(min(i, j), max(i, j)) + 1, gamma);
        }

        double r = double_uniform(engine);

        if (r < p) {
            // update number of interactions
            C(min(i, j), max(i,j)) += 1;

            // choose a random reassignment eps
            double eps = double_uniform(engine);

            // after transaction
            double dm = (1 - lambda)*(eps*m[j] - (1 - eps)*m[i]);
            m[i] += dm;
            m[j] -= dm;
        }

    }

    test_transactions(N, m0, m);


    double mean = sum(m);
    double mean2 = 0;
    for (int i=0; i<N; i++) {
        mean2 += m[i]*m[i];
//        ofile << m[i] << " ";
    }
//    ofile << "\n";

    mean2 /= N;
    mean /= N;

    double var = (mean2 - mean*mean);
    ofile <<txn << " " << var << endl;


}

void test_transactions(int N, double m0, vec &m) {
    double m_init = N*m0;
    double m_actual = sum(m);

    double tol = 1e-5;
    if ((m_init - m_actual) > tol) {
        cout << "Money is not conserved" << endl;
        cout << m_init - m_actual << endl;
    }
}
