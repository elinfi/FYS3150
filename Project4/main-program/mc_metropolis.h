#ifndef MC_METROPOLIS_H
#define MC_METROPOLIS_H

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;


double E_sum (mat A, int L);
int sum_neighbour (mat A, int i, int j, int L);
tuple <double, double, double, double, double> markov_chain (int N, int L, double temp);

#endif // MC_METROPOLIS_H
