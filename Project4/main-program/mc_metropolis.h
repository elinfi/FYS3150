#ifndef MC_METROPOLIS_H
#define MC_METROPOLIS_H

#include <armadillo>

using namespace arma;


double E_sum (mat A, int L);
int sum_neighbour (mat A, int i, int j, int L);
void markov_chain (mat A, int N, int L, double temp);

#endif // MC_METROPOLIS_H
