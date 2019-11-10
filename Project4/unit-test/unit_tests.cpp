#include "test_e_sum.h"
#include "catch.hpp"
#include "../main-program/ising_model.h"
#include "../main-program/mc_metropolis.h"
#include <armadillo>

using namespace std;
using namespace arma;


TEST_CASE("Test E_sum") {
    int L = 3;

    Mat<double> A(L, L, fill::ones);
    for(int i=0; i<L; i++){
        A(i,i) = -1;
    }

    double energy = E_sum(A, L);

    REQUIRE(energy == Approx(6).epsilon(0.000000000001));
}

TEST_CASE("Test sum_neighbour") {
    int L = 3;

    Mat<double> A(L, L, fill::ones);
    A(1, 2) *= -1;

    double dE = sum_neighbour(A, 1, 2, L);

    REQUIRE(dE == Approx(-8).epsilon(0.000000000001));
}

TEST_CASE("Test magnetization") {
    int L = 3;

    Mat<double> A(L, L, fill::ones);
    for(int i=0; i<L; i++){
        A(i,i) = -1;
    }

    double magnet = accu(A);

    REQUIRE(magnet == Approx(3).epsilon(0.000000000001));
}
