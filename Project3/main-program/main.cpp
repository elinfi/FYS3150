#include <iostream>
#include "gauss_legendre.h"
#include "gauss_quadrature.h"
#include "improved_gauss_quadrature.h"
#include <cmath>

using namespace std;


int main()
{
    int n = 10;
//    double lambda = 2.0;

//    double sum = gauss_quadrature(n, lambda);
    double sum = improved_gauss_quadrature(n);

    cout << "N = " << n << endl;
//    cout << "lambda = " << lambda << endl;
    cout << sum << endl;

    double pi = acos(-1);
    double exact = (5*pi*pi)/(16*16);

    cout << abs(exact - sum) << endl;

    return 0;
}
