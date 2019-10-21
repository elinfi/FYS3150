#include <iostream>
#include <armadillo>
#include "mpi/mpi.h"
#define ZERO  1.0E-16
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <time.h>
#define EPS 10.0e-10
#define MAXIT 10


using namespace std;

double gammln(double);
double improved_montecarlo_func(double r1_, double r2_, double thet1, double thet2, double phi1, double phi2);
double montecarlo_two_electric_boogaloo(int n, int seed, double & variance);
double montecarlo_func(double lambda, int n, double &variance);


using namespace std;

double legendre_poly(int n, double x){
    //assumes L0(x)=0 ??
    //L_{j+1} = ((2*j +1)*x*L_j(x) - j*L_{j-1}(x))/(j+1)
    //s = L_{j+1}(x)
    //r = L_j(x)
    //t = L_{j-1}(x)
    double r,s,t;
    r = 0; s = 1;
    //uses recusion relation to generate p1 and p2 (nicked from lecture notes)
    for (int m = 0; m<n; m++){
        t = r; r =s;
        s = (2*m+1)*x*r - m*t;
        s /= (m+1);
    }
    return s;

}

//-----------------------------------a------------------------------------------
void gauleg(double x1, double x2, double x[], double w[], int n){

    int         m;
    double      z1,z,xm,xl,pp,p2,p1;
    double      *x_low, *x_high, *w_low, *w_high;

    m  = (n + 1)/2;                             // roots are symmetric in the interval
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);

    x_low  = x;                                       // pointer initialization
    x_high = x + n - 1;
    w_low  = w;
    w_high = w + n - 1;

    for (int i = 1; i <= m; i++){
        z = cos(M_PI*(i - 0.25)/(n + 0.5));
        do{
            p1 =1.0;
            p2 =0.0;

          /*
          ** loop up recurrence relation to get the
              ** Legendre polynomial evaluated at x
              */

//        for(j = 1; j <= n; j++) {
//           p3 = p2;
//           p2 = p1;
//           p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
//        }


            p1 = legendre_poly(n, z);
            p2 = legendre_poly(n-1, z);
            pp = n*(z*p1-p2)/(z*z-1.0);
            z1 = z;
            z = z1 - p1/pp;
        }while(fabs(z-z1) > 1.0E-10);

        *(x_low ++) = xm - xl*z;
        *(x_high--) = xm + xl*z;
        *w_low = 2.0*xl/((1.0-z*z)*pp*pp);
        *(w_high--) = *(w_low++);
        //Her går koden fra hver sin ende i arrayet, ettersom nullpunkter skal være symetriske
        //om null. Her blir det mindre presisjon om N er et partall!
    }


}

double int_func(double x1, double x2, double y1, double y2, double z1, double z2){
    double alpha = 2;
    double norm = sqrt( ((x1-x2)*(x1-x2))+ ((y1-y2)*(y1-y2)) + ((z1-z2)*(z1-z2)));
    double r1 = sqrt((x1*x1 + y1*y1 + z1*z1));
    double r2 = sqrt((x2*x2 + y2*y2 + z2*z2));
    if(norm < ZERO){
        return 0;
    }
    else {return exp(-2*alpha*(r1 + r2))*(1/norm);}
}


//--------------------------------------------------b------------------------------------
double laguerre_func(double r1_, double r2_, double thet1, double thet2, double phi1, double phi2){


    double r_one_two;
    double cos_betta = cos(thet1)*cos(thet2) + sin(thet1)*sin(thet2)*cos(phi1 - phi2);
    double root = r1_*r1_ + r2_*r2_ - 2*r1_*r2_*cos_betta;
    if (root < EPS){
          return 0.0;
    }

    r_one_two = (sqrt(root));


    //double dr1dr2 = r1_*r1_*exp(-r1_)*r2_*r2_*exp(-r2_);
    return sin(thet1)*sin(thet2)/(1024*r_one_two);
}

//  Function to set up Gauss-Laguerre integration points and weights from
//  the text Numerical Recipes of Teukolsky et al.


//  Note that you need to call it with a given value of alpha,
// called alf here. This comes from x^{alpha} exp(-x)

void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}

// end function gaulag


double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

// end function gammln


//-----------------------------------------------c----------------------------------------
double montecarlo_func(double lambda, int n, double & variance){
    double a = -lambda;
    double b = lambda;
    int seed = 1337;
    mt19937 engine(seed);
    // Generate distribution
    uniform_real_distribution<double> uniform(a, b);

    double variance_sum= 0;

    double sum = 0;
    double res = 0;
    for(int i =0; i<n; i++){

        double x1 = uniform(engine);
        double x2 = uniform(engine);
        double y1 = uniform(engine);
        double y2 = uniform(engine);
        double z1 = uniform(engine);
        double z2 = uniform(engine);


        double value = int_func(x1, x2, y1, y2, z1, z2);
        sum += value;
        variance_sum += value*value;

    }
    double jacobi = pow((b-a), 6);
    res = (sum*jacobi)/n;
    double v1 = variance_sum/n;
    double v2 = sum/n;

    variance = (v1 - v2*v2)*jacobi;
    cout<<variance<<endl;
    return res;
}



//----------------------------------------d----------------------------------------------


double montecarlo_two_electric_boogaloo( int n, int seed, double & variance){

    mt19937 engine(seed);
    uniform_real_distribution<double> uniform(0,1);
    exponential_distribution<double> exponential(1);

    double variance_sum = 0;
    //int seed =
    double sum = 0;
    double res = 0;
    for(int i=0; i<n; i++){
        double r1 = exponential(engine);
        double r2 = exponential(engine);

        double phi1 = uniform(engine)*(2*M_PI);
        double phi2 = uniform(engine)*(2*M_PI);

        double theta1 = uniform(engine)*(M_PI);
        double theta2 = uniform(engine)*(M_PI);

        double value = improved_montecarlo_func(r1, r2, theta1, theta2, phi1, phi2);
        sum += value;
        variance_sum += value*value;


    }


    double jac = pow((M_PI), 2)*pow((2*M_PI), 2);
    res = (sum*jac)/n;
    cout <<res<< " " <<n <<endl;

    double v1 = variance_sum/n;
    double v2 = sum/n;

    variance = (v1 - pow((v2), 2))*jac;




    return res;
}


double improved_montecarlo_func(double r1_, double r2_, double thet1, double thet2, double phi1, double phi2){


    double r_one_two;
    double cos_betta = cos(thet1)*cos(thet2) + sin(thet1)*sin(thet2)*cos(phi1 - phi2);
    double root = r1_*r1_ + r2_*r2_ - 2*r1_*r2_*cos_betta;
    if (root < EPS){
          return 0.0;
    }

    r_one_two = (sqrt(root));


    //double dr1dr2 = r1_*r1_*exp(-r1_)*r2_*r2_*exp(-r2_);
    return (sin(thet1)*sin(thet2)*r1_*r1_*r2_*r2_)/(1024*r_one_two);
}

#undef EPS
#undef MAXIT
#undef ZERO















