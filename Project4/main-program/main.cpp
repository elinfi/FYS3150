#include <iostream>
#include "mc_metropolis.h"
#include <mpi.h>
#include <fstream>

using namespace std;

int main(int args, char* argv[])
{
    MPI_Init(&args, &argv);
    int L = 2;
    int N = 1000000;
    double temp = 1;
    bool random_spin = true;

    tuple <double, double, double, double, double> values = markov_chain(N, L, temp, random_spin);

    double E, E_variance, M, M_variance, M_abs;

//    E = get<0>(values);
//    E_variance = get<1>(values);
//    M = get<2>(values);
//    M_variance = get<3>(values);
//    M_abs = get<4>(values);

    tie(E, E_variance, M, M_variance, M_abs) = values;

    MPI_Finalize();

    double Cv = E_variance/(temp*temp);         // heat capacity
    double Chi = M_variance/(temp*temp);        // susceptibility

    cout << Cv << " " << Chi << endl;

    string file = "4b";
    ofstream ofile;
    ofile.open(file);


    return 0;
}
