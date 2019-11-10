#include "mc_metropolis.h"
#include "ising_model.h"
#include <mpi.h>
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

void markov_chain (mat A, int N, int L, double temp) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int seed = time(0)+world_rank*10;
//    int seed = 1024;

    // generate engine
    mt19937_64 engine(seed);

    // generate integers with uniform distribution
    uniform_int_distribution<int> int_uniform(0, L-1);

    // generate double with uniform distribution
    uniform_real_distribution<double> double_uniform(0, 1);


    int x, y, dE, num;
    double r, k;
    double energy, expected_energy, energy_squared;
    double magnet, expected_magnet, magnet_squared, magnet_abs;
    double Cv, Chi;

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
    num = 0;
    for (int i=0; i<N; i++) {
        // Loop over all spins, pick a random spin each time
        for (int j=0; j<L*L; j++) {
            x = int_uniform(engine);
            y = int_uniform(engine);
//            A(x, y) *= -1;                      // flip a random selected spin in the lattice

            dE = sum_neighbour(A, x, y, L);     // change in energy level

            r = double_uniform(engine);

            if (r <= w[dE + 8]) {
                num += 1;
                A(x, y) *= -1;

                energy += dE;
                magnet += A(x, y)*2;
            }

//            if (r > w[dE + 8]) {                // reject the new state
//                A(x, y) *= -1;
//            }

//            else {                              // accept the new state
//                energy += dE;
//                magnet += 2*A(x, y);
//                num += 1;
//            }
        }


        // Update expectation values
        expected_energy += energy;
        energy_squared += energy*energy;
        expected_magnet += magnet;
        magnet_squared += magnet*magnet;
        magnet_abs += abs(magnet);
    }

    // Normalize average values
    expected_energy /= double(N);
    energy_squared /= double(N);
    expected_magnet /= double(N);
    magnet_squared /= double(N);
    magnet_abs /= double(N);

    Cv = (energy_squared - pow(expected_energy, 2))/(k*temp*temp);        // heat capacity
    Chi = (magnet_squared - pow(expected_magnet, 2))/(k*temp);         // susceptibility

//    double values [] = {expected_energy, energy_squared, expected_magnet, magnet_squared, Cv, Chi};

    cout << num << " " << expected_energy << " " << magnet_abs << " " << expected_magnet << endl;

}
