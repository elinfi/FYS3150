#include <iostream>
#include <armadillo>

using namespace std;


double f(double v1, double v2, double v3, int n) {
    double h = 1/(n + 1);


    double f_i = -(v3 + v1 - 2*v2)/(h*h);

    return f_i;
    // for (int i=1; i<n-1; i++) {
    //    f[i] = -(v[i+1] + v[i-1] - 2*v[i])/h**2
    //}
}

int main(){
    int N = 10;

    double d[N];
    double d_c[N];
    double b[N];
    double b_c[N];

    for (int i=2; i<N-1; i++) {
        // d_c[i] = d[i] - (a[i-1]*c[i-])/d_c[i-1];
        d_c[i] = 2 - 1/d_c[i-1];
        //b_c[i] = b[i] - (b_c[i-1]*a[i-1])/d_c[i-1];
        b_c[i] = b[i] - b_c[i-1]/d_c[i-1];
    }

return 0;
}


