#include <iostream>
#include <armadillo>
#include <fstream>
#include <string>

using namespace std;


double f(double x){
    double f_i = 100*exp(-10*x);

    return f_i;
}


void for_back_alg(double* a, double* b, double* c, double* b_,
                  double* g, double* g_, double* v, int N) {

    // forward substitution
    for (int i=1; i<N; i++) {
        b_[i] = b[i] - (a[i-1]*c[i-1])/b_[i-1];
        g_[i] = g[i] - (g_[i-1]*a[i-1])/b_[i-1];
    }

    // backward substitution
    // first and last point is zero by boundary conditions
    for (int i=N-2; i>0; i--) {
        v[i] = (g_[i] - a[i]*v[i+1])/b_[i];
    }
}



int main(int argc, char *argv[]){

    //object for output file
    ofstream ofile;

    string filename;
    int N;

    if(argc <= 1){
        cout << "Bad usage: " << argv[0] <<
                " read also filename on same line and size of matrix (N)" << endl;
        exit(1);
    }
    else {
        filename = argv[1]; // first command line argument after the name of the program
        N = atoi(argv[2]);
    }


    double h = 1.0/(N + 1);
    double hh = pow(h, 2);


    // creating vectors for matrix A and right side of the equation and initiating all elements to be zero
    double* a = new double[N] ();      // below diagonal
    double* b = new double[N] ();      // diagonal
    double* c = new double[N] ();      // above diagonal
    double* b_ = new double[N] ();     // changed diagonal
    double* g = new double[N] ();      // right side of equation
    double* g_ = new double[N] ();     // changed right side of equation
    double* v = new double[N] ();      // v vector


    // filling up the the array for a, b, c and g
    for (int i=0; i<N; i++){
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
        g[i] = hh*f(i*h);
    }

    // initiating b_
    b_[0] = b[0];


    // solving the differential equation
    for_back_alg(a, b, c, b_, g, g_, v, N);


    // Open file
    ofile.open(filename);

    // write the result to file
    for (int i=0; i<N; i++) {
        if (i == N - 1) {
            ofile << v[i];
        }
        else {
        ofile << v[i] << ",";
        }

    }

    ofile.close();

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] b_;
    delete[] g;
    delete[] g_;
    delete[] v;

return 0;
}


