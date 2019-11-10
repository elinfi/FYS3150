#include <iostream>
#include "mc_metropolis.h"
#include <mpi.h>

using namespace std;

int main(int args, char* argv[])
{
    MPI_Init(&args, &argv);
    int L = 2;
    int N = 1000000;
    double temp = 1;

    tuple <double, double, double, double, double> values = markov_chain(N, L, temp);

    MPI_Finalize();

    return 0;
}
