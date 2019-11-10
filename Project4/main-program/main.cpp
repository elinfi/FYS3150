#include <iostream>
#include "ising_model.h"
#include "mc_metropolis.h"
#include <mpi.h>

using namespace std;

int main(int args, char* argv[])
{
    MPI_Init(&args, &argv);
    int L = 2;
    int N = 1000000;
    double temp = 1;

    ising_model(L, N, temp);

    MPI_Finalize();

    return 0;
}
