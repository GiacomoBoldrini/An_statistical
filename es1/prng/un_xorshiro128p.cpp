#include <iostream>
#include "xorshiro.h"
#include "prng.h"

int main(int argc, char** argv){

    int N = 10;

    if(argc > 1){

        N = std::stoi(argv[1]);
    }

    xorshiro rndxor;

    for(int i = 0; i < N; i++){

        std::cout << rndxor.rand(-10,10) << std::endl;      
    }

    return 0;
}
