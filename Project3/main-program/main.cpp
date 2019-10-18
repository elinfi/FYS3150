#include <iostream>
#include "gauss_legendre.h"
#include "gauss_quadrature.h"
#include "improved_gauss_quadrature.h"
#include "monte_carlo.h"
#include "improved_monte_carlo.h"
#include <cmath>
#include <random>

using namespace std;


int main()
{
    /*
    // Gauss Quadrature
    int n = 25;
    double lambda = 2.7;
    bool timing = true;

    double sum = gauss_quadrature(n, lambda, timing);

    cout << "N = " << n << endl;
    cout << "lambda = " << lambda << endl;
    cout << sum << endl;

    double pi = acos(-1);
    double exact = (5*pi*pi)/(16*16);

    cout << abs(exact - sum) << endl;
    */


    /*
    // Improved Gauss Quadrature

    bool timing = true;
    for (int n=10; n<20 ; n+=5){
        double sum = improved_gauss_quadrature(n, timing);

        cout << "N = " << n << endl;
        cout << sum << endl;

        double pi = acos(-1);
        double exact = (5*pi*pi)/(16*16);

        cout << abs(exact - sum)/exact << endl;
    }
    */


    // Monte Carlo
    int n = 1000000;
    double lambda = 2.7;
    bool timing = true;

    double sum = monte_carlo(n, lambda, timing);

    cout << "N = " << n << endl;
    cout << "lambda = " << lambda << endl;
    cout << "sum: " << sum << endl;

    double pi = acos(-1);
    double exact = (5*pi*pi)/(16*16);

    cout << "relative error: " << abs(exact - sum)/exact << endl;




    /*
    // Improved Monte Carlo

        int n = 100000;
        bool timing = true;
        double sum = improved_monte_carlo(n, timing);

        cout << "N = " << n << endl;
        cout << "sum: " << sum << endl;

        double pi = acos(-1);
        double exact = (5*pi*pi)/(16*16);

        cout << "relative error: " << abs(exact - sum)/exact << endl;
    */

    return 0;
}
