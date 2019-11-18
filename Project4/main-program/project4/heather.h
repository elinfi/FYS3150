#ifndef HEATHER_H
#define HEATHER_H


#include <iostream>
#include <armadillo>
#include "mpi.h"
#include <fstream>


using namespace std;
using namespace arma;

tuple <double, double, double, double, double, double, double> main_func(int L, double T, int n, bool random_spin);
void main_test(bool test);

#endif // HEATHER_H
