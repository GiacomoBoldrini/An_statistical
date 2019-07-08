#include <iostream>
#include "../prng.h"
#include "../prngGen.h"
#include "../xorshiro/xorshiro.h"
#include "bw_xorshiro128p.h"


int main(int argc, char ** argv){

    int N = 100;

    double m0 = 0.;
    double gamma = 2.;

    if(argc > 1){

        N = std::stoi(argv[1]);
    }

    xorshiroGen bw_inv;

    for(int i = 0; i < N; i++){
    
    std::cout << bw_inv.BW_inverse(m0,gamma) << std::endl;

    }

    return 0;
}
