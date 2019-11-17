#include <iostream>
#include "mc_metropolis.h"
#include <fstream>
#include <cmath>
#include <armadillo>
#include <chrono>
using namespace std;

int main()
{

    int N = 50000;
    bool random_spin = false;
    double E, E_variance, M, M_variance, M_abs, M_abs_var;

    double t_start = 2.1;
    double t_end = 2.4;
    double step = 0.005;
    int n = ceil((t_end - t_start)/step) + 1;

    double t_list[n];
    for(int i=0; i<n; i++){
        t_list[i] = t_start + step*i;
    }

    string file;
    ofstream ofile;
    int L = 100;
    file = "4e_L" + to_string(L) + "_N5e5_21_24_005_.txt";
    ofile.open(file);
    ofile << n << endl;
    #pragma omp parallel for
    for (int i=0; i<n; i++) {
        double t = t_list[i];
        tuple <double, double, double, double, double, double> values = markov_chain(N, L, t, random_spin);
        tie(E, E_variance, M, M_variance, M_abs, M_abs_var) = values;

        double Cv = E_variance/(t*t);         // heat capacity
        double Chi = M_variance/t;            // susceptibility
        double Chi_abs = M_abs_var/t;         // susceptibility with absolute magnetization
        cout << t << " " << E << " " << M_abs << " " << M << " " << Cv << " " << Chi << " " << Chi_abs << endl;
        ofile << t << " " << E << " " << M_abs << " " << M << " " << Cv << " " << Chi << " " << Chi_abs << endl;
    }


    /*
    int N = 50000;
    bool random_spin = false;
    double E, E_variance, M, M_variance, M_abs, M_abs_var;

    double t_start = 2.1;
    double t_end = 2.41;
    double step = 0.1;
    int n = ceil((t_end - t_start)/step) + 1;

    double t_list[n];
    for(int i=0; i<n; i++){
        t_list[i] = t_start + step*i;
    }

    int L = 20;
    // start timing
    auto start = chrono::high_resolution_clock::now();


    for (int i=0; i<n; i++) {
        double t = t_list[i];
        tuple <double, double, double, double, double, double> values = markov_chain(N, L, t, random_spin);
        tie(E, E_variance, M, M_variance, M_abs, M_abs_var) = values;

        double Cv = E_variance;         // heat capacity
        double Chi = M_variance;            // susceptibility
        double Chi_abs = M_abs_var/t;         // susceptibility with absolute magnetization
        cout << t << " " << E << " " << M_abs << " " << M << " " << Cv << " " << Chi << " " << Chi_abs << endl;
    }

    // end timing
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = (finish - start);
    cout << "Time: " << elapsed.count() << " s\n" << endl;
    */


    /*
    double temp = 2.4;
    double beta = 1/temp;
    bool random_spin = true;

    double E, E_variance, M, M_variance, M_abs;
    double exact_E, exact_M_abs;

    exact_E = -(8*cosh(8*beta)/(cosh(8*beta) + 3));
    exact_M_abs = (2*exp(8*beta) + 4)/(cosh(8*beta) + 3);

    cout << exact_E/4 << " " << exact_M_abs/4 << endl;

    for (int N=0; N<30000; N++) {
        tuple <double, double, double, double, double> values = markov_chain(N, L, temp, random_spin);
        tie(E, E_variance, M, M_variance, M_abs, M_abs_var) = values;
        double E_err = abs(exact_E - E)/exact_E;
        double M_err = abs(exact_M_abs - M_abs)/exact_M_abs;
        cout << "N: " << N << " Energy: " << E << " E_err: " << E_err <<  " Magnet: " << M_abs << " M_err: " << M_err <<endl;
    }
    */



    /*
    string file = "4c_T1_1_3e5_1_true_L20.txt";
    ofstream ofile;
    ofile.open(file);

    int L = 20;
    double temp = 1;
    bool random_spin = true;
    double E, E_variance, M, M_variance, M_abs, M_abs_var;

    int start = 10;
    int end = 30000;
    int step = 100;


    ofile << (end-start)/step << endl;
    #pragma omp parallel for
    for (int N=start; N<end; N+=step) {
        tuple <double, double, double, double, double, double> values = markov_chain(N, L, temp, random_spin);
        tie(E, E_variance, M, M_variance, M_abs, M_abs_var) = values;
        ofile << N << " " << E << " " << M_abs << " " << M << endl;
        cout << N << " " << E << " " << M_abs << " " << M << endl;
    }
    */


    return 0;
}
