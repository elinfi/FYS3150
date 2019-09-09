#include <iostream>
#include <armadillo>
#include <fstream>
#include <string>
#include "time.h"
#include "general_thomas.h"
#include "optimized_thomas.h"
#include "lu_decomp.h"

using namespace std;


int main(int argc, char *argv[]){
    int N;

    if(argc <= 1){
        cout << "Bad usage: " << argv[0] <<
                " read also size of matrix (N)" << endl;
        exit(1);
    }
    else {
        N = atoi(argv[1]);
    }

    general_thomas(N);
    optimized_thomas(N);
    lu_decomposition(N);

return 0;
}
