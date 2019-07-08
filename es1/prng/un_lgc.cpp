//g++ un_lgc.cpp ../prng.cpp ../prngGen.cpp -o un_lgc
//or ./compile.sh
#include <iostream>
#include <string>
#include "prng.h"


int main(int argc, char** argv){

    int N = 10;

    if(argc > 1){        
        N = std::stoi(argv[1]); //string to integer
    }

    prng gen;

    for(int i = 0; i < N; i++){

        std::cout << gen.rand(-10,10) << std::endl;          
    }

    return 0;
}
