#include "xorshiro.h"
#include "prng.h"
#include <climits>
#include <iostream>



int main(int argc, char** argv){

    xorshiro gen;

    randtype rnd_num = gen.xorshiro128p();

    std::cout.write(reinterpret_cast<char*>(&rnd_num), sizeof rnd_num);

    while(1){

    rnd_num = gen.xorshiro128p();

    std::cout.write(reinterpret_cast<char*>(&rnd_num), sizeof rnd_num);

    }


    return 0;
}