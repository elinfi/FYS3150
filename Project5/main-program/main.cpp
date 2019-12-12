#include <iostream>
#include "transactions.h"
#include <fstream>

using namespace std;

int main()
{
    int N = 500;
    double m0 = 500;
    double lambda = 0.25;
    double alpha = 1.5;
    double gamma = 0.0;
    int txn = 1e4;
    int mc = 1e3;

//    string filename = "var_N" + to_string(N) + "_m0" + to_string(int(m0)) + "_lambda" + to_string(lambda) + "_txn" + to_string(txn) + "_runs" + to_string(runs) + ".txt";
    string filename = "5d_N500_m500_l025_a15_g0_25e4_1e7.txt";
    cout << filename << endl;
    ofstream ofile;
    ofile.open(filename);

    for (int i=25000; i<1e7; i+=50000) {
        transactions(N, m0, lambda, alpha, gamma, i, ofile);
    }
//    for (int i=0; i<mc; i++) {
//        transactions(N, m0, lambda, alpha, gamma, txn, ofile);
//    }
}

/*
 * Beregn variansen for lambda = 0.25, 0.5, 0.9. - finn beste antall transaksjoner
 * Få svar på hva paramteriseringen i 5c betyr. High-end tails.
 * 5d, hva er "perform the same analysis as previously"?
 * 5d extract the tail og sammenlikn med Pareto dsitribution
 * 5d - kjøre for N=500, 1000 og lambda = 0.25, 0.5, 0.9 og alpha = 0.5, 1.0, 1.5,  2.0???
 * 5e wealth distribution?
*/

