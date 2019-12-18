#include <iostream>
#include "transactions.h"
#include <fstream>

using namespace std;

int main()
{
    int N = 1000;
    double m0 = 500;
    double lambda = 0.9;
    double alpha = 2.0;
    double gamma = 0.0;
    int txn = 1e7;
    int mc = 1e4;

    string filename = "5e_N1000_m500_l09_a2_g0_2.txt";
    cout << filename << endl;
    ofstream ofile;
    ofile.open(filename);

    /*
    int t = 5;
    while (t < 1e8) {
            transactions(N, m0, lambda, alpha, gamma, t, true, false, ofile);

            if (t < 10000) {
                t += 1000;
            }
            else if (t > 10000 && t < 5e5) {
                t += 50000;
            }
            else if (t > 5e5 && t < 1e6) {
                t += 70000;
            }
            else if (t > 1e6 && t < 4e6) {
                t += 200000;
            }
            else if (t > 4e6 && t < 1e7) {
                t += 400000;
            }
            else if (t > 1e7 && 5e7) {
                t += 500000;
            }
            else {
                t += 600000;
            }
    }
    */

    /*
    int t = 5;
    while (t < 1e7) {
            transactions(N, m0, lambda, alpha, gamma, t, true, false, ofile);

            if (t < 10000) {
                t += 500;
            }
            else if (t > 10000 && t < 5e5) {
                t += 2000;
            }
            else if (t > 5e5 && t < 1e6) {
                t += 7000;
            }
            else if (t > 1e6 && t < 4e6) {
                t += 25000;
            }
            else {
                t += 50000;
            }
    }
    */

    /*
    for (int i=0; i<mc; i++) {
        transactions(N, m0, lambda, alpha, gamma, txn, false, true, ofile);
    }
    */
}

