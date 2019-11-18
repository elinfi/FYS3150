#include <iostream>
#include <armadillo>
#include "mpi.h"
// #include <mpi/mpi.h>
#define ZERO  1.0E-10
#include "heather.h"
#include <fstream>
#include <time.h>
#include <random>

using namespace std;
using namespace arma;

void main_test(bool test);


double E_sum(mat A, int L){
    //Sum of nearest neighbours with periodic boundry condition
    int sum = 0;
    for(int i=0; i<L; i++){
        for(int j =0; j<L; j++){
            sum-= A(i,j)*A((i+1)%L, j) +A(i, j)*A(i, (j+1)%L);
        }
    }
    return sum;
}

int sum_neighbour (mat A, int i, int j, int L) {
    int dE = 2*A(i, j)*(A((i+1)%L, j) + A(i, (j+1)%L)+ A((i+L-1)%L, j) + A(i, (j+L-1)%L));

    return dE;
}

tuple <double, double, double, double, double, double, double> main_func(int L, double T, int n, bool random_spin ){
    mat A(L,L, fill::ones );
//    int world_rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int seed = time(0);//+world_rank*10;
    mt19937 engine(seed);

    uniform_real_distribution<double> zero_one(0, 1); //used for finding 'r' the variable we check in Metropolis
    uniform_int_distribution<int> uniform(0, L-1); //used for finding indices in the matrix
    uniform_int_distribution<int> uniform2(0,1); //used in generating matrix


//  filling matrix with random ones or zeros using mercenne twister

    if (random_spin) {
        for (int i=0; i<L; i++) {
            for (int j=0; j<L; j++) {
                int num = uniform2(engine);

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
    int x, y, dE;
    double r, k;
    double energy, expected_energy, energy_squared, E_variance;
    double magnet, expected_magnet, magnet_squared, M_variance, magnet_abs;
    tuple <double, double, double, double, double, double, double> values;

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
        w[de + 8] = exp(-de/T);
    }


    dE = 0;

    int acceptance = 0;
    ofstream myfile;
//    myfile.open(to_string(L)+"-"+to_string(T)+"-"+to_string(random_spin)+"-"+"ENERGIES_acceptance"+".txt");
    for (int i=0; i<n; i++) {
        // Loop over all spins, pick a random spin each time
        for (int j=0; j<L*L; j++) {

            x = uniform(engine);
            y = uniform(engine);
//            A(x, y) *= -1;                      // flip a random selected spin in the lattice

            dE = sum_neighbour(A, x, y, L);     // change in energy level

            r = zero_one(engine);


            if (r <= w[dE + 8]) {
                A(x, y) *= -1.0;

                energy += dE;
                acceptance ++;

                magnet += A(x, y)*2;
            }//end if r<= w




        }
//        if(i > 7000){
//            myfile << to_string(energy)+ " " + to_string(acceptance)+"\n"; //writing energies to file for plotting of histogram
            // Update expectation values
            expected_energy += energy;
            energy_squared += energy*energy;
            expected_magnet += magnet;
            magnet_squared += magnet*magnet;
            magnet_abs += abs(magnet);
//        }//end if i<7000
    }//end for i<n

    // Normalize average values
    expected_energy /= double(n);
    energy_squared /= double(n);
    expected_magnet /= double(n);
    magnet_squared /= double(n);
    magnet_abs /= double(n);

    E_variance = (energy_squared - pow(expected_energy, 2))/(L*L*T*T);  //if printing heat capacity, don't dive by T.
    M_variance = (magnet_squared - pow(expected_magnet, 2))/(L*L*T); // dividing per spins

    expected_energy /= double(L*L); //per spin
    expected_magnet /= double(L*L);
    magnet_abs /= double(L*L);
    double M_abs_variance = (magnet_squared - pow(magnet_abs, 2)/(L*L*T));


    values = make_tuple(expected_energy, E_variance, expected_magnet, M_variance, magnet_abs, acceptance, M_abs_variance);





    return values;
}//end function











void main_test(bool test){
    if(test){



        double E_exp = 0; double E_variance;  double M_exp = 0;double M_variance = 0;
        double acceptance = 0; double M_abs = 0; double M_abs_variance = 0;

        int L = 2;
        int n=50000000;

        double tol = 1E-3;
        double anal = -1.99598;
        //part 4C
        double T = 1;


        //benchmarc call for T=1, Kb = 1, L = 2
        tuple <double, double, double, double, double, double, double> values = main_func(L, T, n, true);

        tie(E_exp, E_variance, M_exp, M_variance, M_abs, acceptance, M_abs_variance) = values;

        cout << E_exp << M_exp<< endl;
        if(fabs(E_exp - anal ) > tol){
            cout<<E_exp<<endl;
            cout << anal << endl;
            cout << E_exp - anal <<endl;
            cout << "Error in unit test \n";
            abort();
        }



    }//end if(test)

}//end main_test
