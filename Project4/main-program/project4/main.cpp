#include <iostream>
#include <armadillo>
#include "mpi.h"
//#include <mpi/mpi.h>
#define ZERO  1.0E-10
#include "heather.h"
#include <fstream>
#include <time.h>
#include <random>


using namespace std;
using namespace arma;
void c_call();
void d_call();
void e_call();

int main(int argc, char* argv[])
{

//    main_test(true);
//    MPI_Init(&argc, &argv);

    c_call();
    //d_call();
    //e_call();





}

void e_call(){
    double E_exp = 0; double E_variance;  double M_exp = 0;double M_variance = 0;
    double acceptance = 0; double M_abs = 0; double M_abs_variance = 0;

    int L = 2;
    int n = 50000;
    bool random_spins = false;
    ofstream myfile;
    float T_start = 2.0;
    float T_end = 2.5;

    string s0 = "-";
    string s1 = to_string(T_start);
    string s2 = to_string(T_end);
    string s3 = to_string(random_spins);

    myfile.open(s1+s0+s2+s0+s3+"e"".txt");

    for(float T = 2.1; T < 2.4; T+=0.1){
        tuple <double, double, double, double, double, double, double> values = main_func(L, T, n, random_spins);

        tie(E_exp, E_variance, M_exp, M_variance, M_abs,  acceptance, M_abs_variance) = values;
        myfile << to_string(T) << " " << to_string(E_exp) << " " <<to_string(E_variance) << " " << to_string(M_exp) << " " << to_string(M_variance) << " " << to_string(M_abs)<< " " << to_string(acceptance) << " " << to_string(M_abs_variance) <<"\n";

    }//end T loop
    myfile.close();

}//end e_call


void d_call(){
    double E_exp = 0; double E_variance;  double M_exp = 0;double M_variance = 0;
    double acceptance = 0; double M_abs = 0; double M_abs_variance = 0;
    ofstream myfile;

    int L = 20;
    int n = 200000;
    bool random_spins = true;
    double T = 1;

    string s1 = to_string(L);
    string s2 = to_string(T);
    string s3 = to_string(random_spins);
    cout<<s3<<endl;
    tuple <double, double, double, double, double, double, double> values = main_func(L, T, n, random_spins);
    tie(E_exp, E_variance, M_exp, M_variance, M_abs, acceptance, M_abs_variance) = values;

    myfile.open(s1 +"-"+ s2+"-"+s3+".txt");


    myfile << to_string(E_exp) << " " << to_string(M_variance) << " " << to_string(M_exp) << " " << to_string(M_variance) << " " <<to_string(acceptance)<< " " <<to_string(n)<<"\n";




}

void c_call(){
    double E_exp = 0; double E_variance;  double M_exp = 0;double M_variance = 0;
    double acceptance = 0; double M_abs = 0; double M_abs_variance = 0;
    bool random_spin = false;


    ofstream myfile;
    int L = 20;
    string s1 = to_string(L);


    int n = 8000;
    //part 4C



    for(float T = 1; T<2.5; T+=1.4){
        for(int j = 0; j<2; j++){


            string s2 = to_string(T);
            string s3 = to_string(random_spin);
            string s4 = to_string(n);
            string s0 = " ";

            myfile.open(s1+s0 +s2+ s0 + s3 +s0 + s4 + s0 + "magnet_plot"+".txt");

            //benchmarc call for T=1, Kb = 1, L = 2
            tuple <double, double, double, double, double, double, double> values = main_func(L, T, n, random_spin);
    //        MPI_Finalize();
            tie(E_exp, E_variance, M_exp, M_variance, M_abs, acceptance, M_abs_variance) = values;


            cout << E_variance << "  " << M_variance<<endl;
            myfile << to_string(E_exp) << " " << to_string(E_variance) << " " << to_string(M_exp) << " " << to_string(M_variance) << " " <<to_string(M_abs)<< " " <<to_string(n)<<"\n";


            random_spin = true;
        }
        }
    myfile.close();
}


