#include <iostream>
#include "gauss_legendre.h"
#include "gauss_quadrature.h"
#include "improved_gauss_quadrature.h"
#include "monte_carlo.h"
#include "improved_monte_carlo.h"
#include <cmath>
#include <random>
#include <mpi.h>
#include <chrono>
#include <tuple>
#include <fstream>

using namespace std;

tuple<double, double> print_gauss_quadrature(int n, double lambda, bool timing) {
    double sum, exact, rel_err;
    cout << "Gauss Quadrature" << endl;

    cout << "N = " << n << endl;
    cout << "lambda = " << lambda << endl;

    // start timing
    auto start = chrono::high_resolution_clock::now();

    sum = gauss_quadrature(n, lambda);

    // end timing
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = (finish - start);

    cout << "Sum: " << sum << endl;

    double pi = acos(-1);
    exact = (5*pi*pi)/(16*16);
    rel_err = abs(exact - sum)/exact;

    cout << "Relative error: " << rel_err << endl;

    if (timing) {
        // print time
        cout << "Time: " << elapsed.count() << " s\n" << endl;
    }

    return make_tuple(rel_err, elapsed.count());
}

tuple<double, double> print_improved_gauss_quadrature(int n, bool timing) {
    double sum, exact, rel_err;

    cout << "Improved Gauss Quadrature" << endl;

    // start timing
    auto start = chrono::high_resolution_clock::now();

    sum = improved_gauss_quadrature(n);

    // end timing
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = (finish - start);

    cout << "N = " << n << endl;
    cout << "Sum: " << sum << endl;

    double pi = acos(-1);
    exact = (5*pi*pi)/(16*16);
    rel_err = abs(exact - sum)/exact;

    cout << "Relative error: " << rel_err << endl;

    if (timing) {
        // print time
        cout << "Time: " << elapsed.count() << " s\n" << endl;
    }

    return make_tuple(rel_err, elapsed.count());
}

tuple<double, double> print_monte_carlo(int n, double lambda, bool timing, int seed) {
    double sum, variance, exact, rel_err;

    cout << "Brute force Monte Carlo" << endl;

    cout << "N = " << n << endl;
    cout << "lambda = " << lambda << endl;

    // start timing
    auto start = chrono::high_resolution_clock::now();

    tie(sum, variance) = monte_carlo(n, lambda, seed);

    // end timing
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = (finish - start);

    cout << "sum: " << sum << endl;

    double pi = acos(-1);
    exact = (5*pi*pi)/(16*16);
    rel_err = abs(exact - sum)/exact;

    cout << "relative error: " << rel_err << endl;
    cout << "variance: " << variance << endl;

    if (timing) {
        // print time
        cout << "Time: " << elapsed.count() << " s\n" << endl;
    }

    return make_tuple(rel_err, elapsed.count());
}

tuple<double, double> print_improved_monte_carlo(int n, bool timing, int seed) {
    double sum, variance, exact, rel_err;

    cout << "Improved Monte Carlo" << endl;

    // start timing
    auto start = chrono::high_resolution_clock::now();

    tie(sum, variance) = improved_monte_carlo(n, seed);

    // end timing
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = (finish - start);

    cout << "N = " << n << endl;
    cout << "sum: " << sum << endl;

    double pi = acos(-1);
    exact = (5*pi*pi)/(16*16);
    rel_err = abs(exact - sum)/exact;

    cout << "relative error: " << rel_err << endl;

    if (timing) {
        // print time
        cout << "Time: " << elapsed.count() << " s\n" << endl;
    }

    return make_tuple(rel_err, elapsed.count());
}

tuple<double, double> print_parallellized_monte_carlo(int &nargs, char* args[], int n, double lambda, bool timing) {
    int numprocs, my_rank, flag;
    double time_start, time_end, total_time;
    double total_sum, local_sum, variance, rel_err, exact;

    //   MPI initializations
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // start timing
    time_start = MPI_Wtime();

    total_sum = 0.0;
    tie(local_sum, variance) = monte_carlo(n, lambda, time_t(0) + my_rank);

    MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    total_sum = total_sum/numprocs;

    // end timing
    time_end = MPI_Wtime();
    total_time = time_end-time_start;

    if (my_rank == 0) {
        cout << "Parallellized Monte Carlo" << endl;
        cout << "N = " << n << endl;
        cout << "Lambda = " << lambda << endl;
        cout << "Sum: " << total_sum << endl;

        double pi = acos(-1);
        exact = (5*pi*pi)/(16*16);
        rel_err = abs(exact - total_sum)/exact;

        cout << "Relative Error: " << rel_err << endl;
        cout << "Variance: " << variance << endl;

        if (timing) {
            cout << "Timing: " << total_time << " s\n" << endl;
        }
    }

    //  End MPI
//    MPI_Finalize ();

    return make_tuple(rel_err, total_time);
}

tuple<double, double> print_parallellized_improved_monte_carlo(int &nargs, char* args[], int n, bool timing) {
    int numprocs, my_rank, flag;
    double time_start, time_end, total_time;
    double total_sum, local_sum, variance, rel_err, exact;

    //   MPI initializations
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // start timing
    time_start = MPI_Wtime();

    // Parallellized Monte Carlo
    total_sum = 0.0;
    tie(local_sum, variance) = improved_monte_carlo(n, time_t(0) + my_rank);

    MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    total_sum = total_sum/numprocs;

    // end timing
    time_end = MPI_Wtime();
    total_time = time_end-time_start;

    if (my_rank == 0) {
        cout << "Parallellized Monte Carlo" << endl;
        cout << "N = " << 2*n << endl;
        cout << "Sum: " << total_sum << endl;

        double pi = acos(-1);
        exact = (5*pi*pi)/(16*16);
        rel_err = abs(exact - total_sum)/exact;

        cout << "Relative Error: " << rel_err << endl;
        cout << "Variance: " << variance << endl;

        if (timing) {
            cout << "Timing: " << total_time << " s\n" << endl;
        }
    }

    //  End MPI
//    MPI_Finalize ();

    return make_tuple(rel_err, total_time);
}


int main(int nargs, char* args[]) {
    int n;
    double lambda;
    bool timing;

    lambda = 1.9;
    timing = true;

    double rel_err, time, rel_err_improved, time_improved;

    string filename = "parallellized_improved_mc_O0";
    ofstream ofile;
    ofile.open(filename);

    /*
    for (n = 1; n < 40; n += 2) {
        tie(rel_err, time) = print_gauss_quadrature(n, lambda, timing);
        ofile << rel_err << " " << time << "\n";
    }
    */


    /*
    n = 100;
    while (n < 100000000) {
        tie(rel_err, time) = print_parallellized_improved_monte_carlo(nargs, args, n, timing);
        ofile << rel_err << " " << time << "\n";

        if (n < 100000) {
            n += 2000;
        }
        else if (n < 1000000) {
            n += 200000;
        }
        else if (n < 10000000){
            n += 2000000;
        }
        else {
            n += 20000000;
        }
    }
    */


    /*
    for (n=100; n < 100000000; n *= 10) {
        tie(rel_err, time) = print_parallellized_improved_monte_carlo(nargs, args, n, timing);
        ofile << 2*n << " " << time << "\n";
    }

    ofile.close();
    MPI_Finalize ();
    */


    /*
    for (n=5; n<40; n += 5) {
        tie(rel_err, time) = print_gauss_quadrature(n, lambda, timing);
        tie(rel_err, time) = print_improved_gauss_quadrature(n, timing);
    }
    */


    /*
    for (n=100; n < 10000000000; n *= 10) {
        tie(rel_err, time) = print_improved_monte_carlo(n, timing, time_t(0));
    }
    */

    return 0;
}
