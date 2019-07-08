#include <iostream>
#include "prng.h"
#include "prngGen.h"


int main(int argc, char ** argv){

    int N = 100;

    double m0 = 0.;
    double gamma = 2.;

    if(argc > 1){

        N = std::stoi(argv[1]);
    }

    prngGen bw_hm;

    for(int i = 0; i < N; i++){
    
    std::cout << bw_hm.BW_hm(m0,gamma) << std::endl;

    }

    return 0;
} 
