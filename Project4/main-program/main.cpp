#include <iostream>
#include "mc_metropolis.h"
//#include <mpi.h>
#include <fstream>
#include <cmath>
#include <armadillo>
using namespace std;

int main()
{
//    int numprocs, my_rank;

//    MPI_Init(&args, &argv);
//    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
//    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    int N = 50000;
//    int L = 40;
    bool random_spin = true;
    double E, E_variance, M, M_variance, M_abs;

    double t_start = 2.2;
    double t_end = 2.4;
    double step = 0.01;
    int n = ceil((t_end - t_start)/step) + 1;

    double t_list[n];
    for(int i=0; i<n; i++){
        t_list[i] = t_start + step*i;
        cout << t_list[i] << endl;
    }

    string file;
    ofstream ofile;
    ofile.open(file);
    ofile << n << endl;

    #pragma omp parallel for
    for (int L = 40; L < 110; L+=20) {
        file = "4e_" + to_string(L) + "_N5e5_22_24_01.txt";
        ofile.open(file);
        ofile << n << endl;
        for (int i=0; i<n; i++) {
            double t = t_list[i];
            tuple <double, double, double, double, double> values = markov_chain(N, L, t, random_spin);
            tie(E, E_variance, M, M_variance, M_abs) = values;

            double Cv = E_variance/(t*t);         // heat capacity
            double Chi = M_variance/(t*t);        // susceptibility
            cout << t << " " << E << " " << M_abs << " " << Cv << " " << Chi << endl;
            ofile << t << " " << E << " " << M_abs << " " << Cv << " " << Chi << endl;
        }
    }




//    double temp = 2.4;
//    double beta = 1/temp;
//    bool random_spin = true;

//    double E, E_variance, M, M_variance, M_abs;
//    double exact_E, exact_M_abs;

//    exact_E = -(8*cosh(8*beta)/(cosh(8*beta) + 3));
//    exact_M_abs = (2*exp(8*beta) + 4)/(cosh(8*beta) + 3);

//    cout << exact_E/4 << " " << exact_M_abs/4 << endl;

    /*
    for (int N=100; N<100000; N*=10) {
        tuple <double, double, double, double, double> values = markov_chain(N, L, temp, random_spin);
        tie(E, E_variance, M, M_variance, M_abs) = values;
        double E_err = abs(exact_E - E)/exact_E;
        double M_err = abs(exact_M_abs - M_abs)/exact_M_abs;
        cout << "N: " << N << " Energy: " << E << " E_err: " << E_err <<  " Magnet: " << M_abs << " M_err: " << M_err <<endl;
    }
    */




//    cout << Cv << " " << Chi << endl;





    /*
    string file = "4c_T24_10_5e5_300_true.txt";
    ofstream ofile;
    ofile.open(file);

    int start = 10;
    int end = 50000;
    int step = 300;
    ofile << (end-start)/step << endl;
    for (int N=start; N<end; N+=step) {
        tuple <double, double, double, double, double> values = markov_chain(N, L, temp, random_spin);
        tie(E, E_variance, M, M_variance, M_abs) = values;
        ofile << N << " " << E << " " << M_abs << endl;
    }
    */

//    MPI_Finalize();


    return 0;
}
