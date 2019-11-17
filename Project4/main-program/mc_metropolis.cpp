#include "mc_metropolis.h"
#include "visualize.h"
#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;


double E_sum (mat A, int L) {
    int sum = 0;
    for (int i=0; i<L; i++) {
        for (int j=0; j<L; j++) {
            sum -= A(i, j)*(A((i+1)%L, j) + A(i, (j+1)%L));
        }
    }
    return sum;
}

int sum_neighbour (mat A, int i, int j, int L) {
    int dE = 2*A(i, j)*(A((i+1)%L, j) + A(i, (j+1)%L)+ A((i+L-1)%L, j) + A(i, (j+L-1)%L));

    return dE;
}

tuple <double, double, double, double, double, double> markov_chain (int N, int L, double temp, bool random_spin) {
    int seed = time(0) + omp_get_thread_num()*10;

    // generate engine
    mt19937_64 engine(seed);

    // generate integers in interval [0, 1] with uniform distribution
    uniform_int_distribution<int> uniform(0, 1);

    // generate integers in interval [0, L-1] with uniform distribution
    uniform_int_distribution<int> int_uniform(0, L-1);

    // generate double with uniform distribution
    uniform_real_distribution<double> double_uniform(0, 1);


    int x, y, dE;
    double r, k;
    double energy, expected_energy, energy_squared, E_variance;
    double magnet, expected_magnet, magnet_squared, M_variance, magnet_abs, M_abs_variance;
    tuple <double, double, double, double, double, double> values;

    // Initializing spin matrix
    Mat<double> A(L, L, fill::ones);

    if (random_spin) {
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
    }

    // initialize variables
    k = 1;                  // Boltzmann constant

    // energy at first microstate
    energy = E_sum(A, L);

    // magnetic moment at first microstate
    magnet = accu(A);       // sum up all spins

    // energy
    expected_energy = energy;
    energy_squared = pow(expected_energy, 2);

    // magnetic moment
    expected_magnet = magnet;
    magnet_abs = abs(magnet);
    magnet_squared = magnet*magnet;

    // vector with the possible energy states
    double* w = new double[17];
    for (int de=-8; de < 9; de+=4) {
        w[de + 8] = exp(-de/temp);
    }

    dE = 0;
    for (int i=0; i<N; i++) {
        // Loop over all spins, pick a random spin each time
        for (int j=0; j<L*L; j++) {
            x = int_uniform(engine);
            y = int_uniform(engine);

            dE = sum_neighbour(A, x, y, L);     // change in energy level

            r = double_uniform(engine);

            if (r <= w[dE + 8]) {
                A(x, y) *= -1.0;

                energy += dE;
                magnet += A(x, y)*2;
            }
        }


        if (N > 7000) {
            // Update expectation values
            expected_energy += energy;
            energy_squared += energy*energy;
            expected_magnet += magnet;
            magnet_squared += magnet*magnet;
            magnet_abs += abs(magnet);
        }
    }

    // Normalize average values
    expected_energy /= double(N);
    energy_squared /= double(N);
    expected_magnet /= double(N);
    magnet_squared /= double(N);
    magnet_abs /= double(N);

    E_variance = (energy_squared - pow(expected_energy, 2))/(L*L);
    M_variance = (magnet_squared - pow(expected_magnet, 2))/(L*L);
    M_abs_variance = (magnet_squared - pow(magnet_abs, 2))/(L*L);

    expected_energy /= double(L*L);
    expected_magnet /= double(L*L);
    magnet_abs /= double(L*L);

    values = make_tuple(expected_energy, E_variance, expected_magnet, M_variance, magnet_abs, M_abs_variance);

    return values;
}
