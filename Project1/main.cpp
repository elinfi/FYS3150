#include <iostream>
#include <armadillo>
#include <fstream>
#include <string>

using namespace std;



//double f(double v1, double v2, double v3, int n) {

// x_1 = i*h
// f(x) = 100*exp(-10*x)
// b = pow(h, 2)*f_i


double f(double x){
    double f_i = 100*exp(-10*x);

    return f_i;
}


//double f(double *v[], int N){

//    double h = 1/(N + 1);
//    double f[N-2];
//    double h_2 = 1/(pow(h, 2));

//    double f_i = -(v3 + v1 - 2*v2)/(pow(h, 2));

//    for (int i=1; i<N-1; i++) {
//        f[i-1] = -(v[i+1] + v[i-1] - v[i])*h_2;


int main(int argc, char *argv[]){
    //object for output file
    ofstream ofile;

    string filename;
    int N;

    if(argc <= 1){
        cout << "Bad usage: " << argv[0] <<
                " read also filename on same line and size of matrix" << endl;
        exit(1);
    }
    else {
        filename = argv[1]; // first command line argument after the name of the program
        N = atoi(argv[2]);
    }


    double h = 1.0/(N + 1);
    double hh = pow(h, 2);

    double d_c[N];
    double b[N];
    double b_c[N];
    double u[N];


    for (int i=0; i<N; i++){
        b[i] = hh*f(i*h);
        cout << b[i] << endl;
    }


    d_c[0] = 2;
    d_c[N] = 0;
    b_c[0] = b[0];


//    d = new double[N];
//    d_c = new double[N];
//    b = new double[N];
//    b_c = new double[N];




    for (int i=2; i<N-1; i++) {
        // d_c[i] = d[i] - (a[i-1]*c[i-])/d_c[i-1];
        d_c[i] = 2 - 1/d_c[i-1];
        //b_c[i] = b[i] - (b_c[i-1]*a[i-1])/d_c[i-1];
        b_c[i] = b[i] + b_c[i-1]/d_c[i-1];
    }

    for (int j=N-2; j>0; j--) {
        u[j] = (b_c[j] + u[j+1])/d_c[j];
    }


    // Declare new filename
    string fileout = filename;
    // Convert the size N to string
    string argument = to_string(N);
    // Final filename as filename-N
//    fileout.append(argument);
    cout << fileout << endl;
    // Open file
    ofile.open(fileout);


    for (int i=0; i<N; i++) {
//        ofile << solution[i] << ",";

        if (i == N - 1) {
            ofile << u[i];
        }
        else {
        ofile << u[i] << ",";
        }

    }

    ofile.close();

//    delete[] d;
//    delete[] d_c;
//    delete[] b;
//    delete[] b_c;

return 0;
}


