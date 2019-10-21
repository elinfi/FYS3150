#include <iostream>
#include <armadillo>
#include <mpi.h>
#include <mpi/mpi.h>
#define ZERO  1.0E-10
#include "header.h"
#include <fstream>
#include <time.h>

using namespace std;
//using namespace arma;

double gauleg_call(int n, double lambda);
void lagrange_call_call();
void gauleg_call_call();
double lagrange_call(int n);
void simple_montecarlo_call();
double rel_err(double numeric);
void montecarlo_two_electric_boogaloo_call();


int main(int argc, char* argv[]) {

    time_t start, stop;


    /*
    //Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//    cout << world_rank << endl;

    */

    start = clock();



    //montecarlo_func(1.9, 1000000);
    stop = clock();


    //lagrange_call(20);
    montecarlo_two_electric_boogaloo_call();
//    MPI_Finalize();

    double time_used = (double) (stop - start)/(CLOCKS_PER_SEC);
    //cout << time_used << "s" <<endl;

    //lagrange_call_call();
    //simple_montecarlo_call();

    return 0;
}

//------------------------------------a-------------------------------------------------

void gauleg_call_call(){
    time_t start, stop;
    string s1, s2;
    double res;
    ofstream myfile;
    myfile.open("a.txt");
    int n;
    int m = 50;
    double lambda;


    for (lambda = 1.0; lambda < 3; lambda +=0.3){
        s2 = to_string(lambda);
        for(n = 5; n < m; n+= 7){
            start = clock();
            res = gauleg_call(n, lambda);
            stop = clock();
            double time_used = (double) (stop - start)/(CLOCKS_PER_SEC);
            double err = rel_err(res);
            myfile << to_string(res) << " " << to_string(n) << " " << to_string(lambda)<< " " << to_string(time_used) << " " << to_string(err) << "\n";

        }
        myfile << "--------------------------------- \n";
    }
    myfile.close();


}

double gauleg_call(int n, double lambda){
    double sum = 0;
    double omega;
    double *x = new double[n];
    double *w = new double[n];

    gauleg( -lambda, lambda, x, w, n);


    for(int i = 0; i<n;i++){
        for(int j = 0; j <n; j++){
            for(int k=0; k < n; k++){
                for(int l=0; l<n; l++){
                    for(int m = 0; m<n; m++){
                        for(int o = 0; o<n; o++){
                            omega = w[i]*w[j]*w[k]*w[l]*w[m]*w[o];
                            sum += omega*int_func(x[i], x[j], x[k], x[l], x[m], x[o]);

                        }
                    }
                }
            }
        }
    }
    cout << sum<<endl;




    delete [] x;
    delete [] w;
    return sum;

}


//-----------------------------------------b-------------------------------------------

void lagrange_call_call(){
    time_t start, stop;
    string s1, s2;
    double res;
    ofstream myfile;
    myfile.open("b.txt");
    int n;
    int m = 40;



    for(n = 5; n < m; n+= 5){
        start = clock();
        res = lagrange_call(n);
        stop = clock();
        double time_used = (double) (stop - start)/(CLOCKS_PER_SEC);
        double err = rel_err(res);
        myfile << to_string(res) << " " << to_string(n) << " " << to_string(time_used) << " " << to_string(err) << "\n";


    }
    myfile.close();


}



double lagrange_call(int n){
    double sum = 0;
    double omega;
    double *x = new double[n+1];
    double *w = new double[n+1];
    double alf = 2;

    double *w_phi = new double[n];
    double *x_phi = new double[n];

    double *w_theta = new double[n];
    double *x_theta = new double[n];

    gauleg(0, M_PI, x_theta, w_theta, n);
    gauleg(0, 2*M_PI, x_phi, w_phi, n);

    gauss_laguerre(x, w, n, alf);



    for(int i = 1; i<(n+1);i++){ //phi1
        for(int j = 1; j <(n+1); j++){ //theta1
            for(int k=0; k < n; k++){ //phi2
                for(int l=0; l<n; l++){ //theta2
                    for(int m = 0; m<n; m++){ //r1
                        for(int o = 0; o<n; o++){ //r2
                            omega = w[i]*w[j]*w_theta[k]*w_theta[l]*w_phi[m]*w_phi[o];
                            sum += omega*laguerre_func(x[i], x[j], x_theta[k], x_theta[l], x_phi[m], x_phi[o]);

                        }
                    }
                }
            }
        }
    }
    cout << sum<<endl;
    delete [] x;
    delete [] w;
    delete [] x_phi;
    delete [] w_phi;
    delete [] x_theta;
    delete [] w_theta;
    return sum;
}

//-----------------------------------------------c-------------------------------------

void simple_montecarlo_call(){
    double res;

    ofstream myfile;

    myfile.open("c.txt");
    int start = 10;
    int stop = 50000;
    int step = 1000;
    double variance =0;
    string s1, s2, s3, s4;
    int n = ((int)stop - start)/(step -1) +1;
    cout << n << endl;
    s4 = to_string(n);
    myfile << s4 + " \n";
    for(double lambd = 1.3; lambd < 4; lambd += 0.3){
        s2 = to_string(lambd);
        for(int i=10; i<stop; i += step){
            s1 = to_string(i);
            res = montecarlo_func(lambd, i, variance);
            s3 = to_string(res);
            s4 = to_string(variance);
            myfile << s3+" "+s2+" "+s1+" " + s4 + "\n";

        }
        myfile <<"------------------ \n";

    }
    myfile.close();


    cout<<res<<endl;
}

double rel_err(double numeric){
    double anal = (5*M_PI*M_PI)/(16*16);
    return((numeric - anal)/anal);

}

void montecarlo_two_electric_boogaloo_call(){
    double res;

    ofstream myfile;

    myfile.open("boogaloo.txt");
    int start = 10;
    int stop = 1E7;
    int step = 20000;
    double mont_var =0;
    int mont_seed = 666;
    string s1, s2, s3, s4;
    int n = ((int)stop - start)/(step -1) +1;
    cout << n << endl;
    s4 = to_string(n);
    myfile << s4 + " \n";
    for(int i=10; i<stop; i += step){
        s1 = to_string(i);
        res = montecarlo_two_electric_boogaloo(i, mont_seed, mont_var);
        s3 = to_string(res);
        s4 = to_string(mont_var);
        myfile << s3+" "+s2+" "+s1+" " + s4 + "\n";

    }
    myfile <<"------------------ \n";


    myfile.close();


    cout<<res<<endl;



}
#undef ZERO
