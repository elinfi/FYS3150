#include <iostream>
#include <math.h>

double analytic_solution (double x) {
        return (1.0 - (1- exp(-10))*x - exp(-10*x));
}

double f (double x){
    return (100*exp(-10*x));
}
