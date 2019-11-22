#include <iostream>
#include <transactions.h>
#include <fstream>

using namespace std;

int main()
{
    int N = 500;
    double m0 = 500;
    int txn = 1e7;
    int runs = 10;

    string filename = "5a";
    ofstream ofile;
    ofile.open(filename);

    for (int i=0; i<runs; i++) {
        transactions(N, m0, txn, ofile);
    }
}
