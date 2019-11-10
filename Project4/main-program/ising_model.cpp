#include "ising_model.h"
#include "mc_metropolis.h"
#include <armadillo>
#include <mpi.h>
#include <iostream>
#include <random>

using namespace std;
using namespace arma;


void ising_model (int L, int N, double temp) {
    Mat<double> A(L, L, fill::zeros);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int seed = time(0)+world_rank*10;
//    int seed = 1024;

    // generate engine
    mt19937_64 engine(seed);

    // generate uniform distribution
    uniform_int_distribution<int> uniform(0, 1);

    for (int i=0; i<L; i++) {
        for (int j=0; j<L; j++) {
            int num = uniform(engine);

            if (num == 0) {
                A(i, j) = -1;
            }
            else {
                A(i, j) = num;
            }
        }
    }

    markov_chain(A, N, L, temp);
}

