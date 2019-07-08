#include <iostream>
#include "prng.h"
#include "prngGen.h"
#include "xorshiro.h"
#include "xorshiroGen.h"



int main(int argc, char ** argv){

    int N = 100;

    double m0 = 2.;
    double gamma = 0.5;

    if(argc > 1){

        N = std::stoi(argv[1]);
    }

    xorshiroGen landau;

    for(int i = 0; i < N; i++){
    
    std::cout << landau.landau_tc(m0,gamma,0,10,0,10) << std::endl;

    //std::cout << landau.l(m0,gamma,0,10,0,10) << std::endl;

    }

    return 0;
}
