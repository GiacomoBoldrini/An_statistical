#include "xorshiro.h"
#include "prng.h"
#include <cmath>
#include <climits>


xorshiro::xorshiro(){

    a = 53;
    b = 22;
    c = 41;

    s[0] = rnd.randInt();
    s[1] = rnd.randInt();

   }

xorshiro::~xorshiro(){}


randtype xorshiro::rotl(const randtype x, int k){

    return (x << k) | (x >> (64 - k ));

}

randtype xorshiro::xorshiro128p(){

    const randtype s0 = s[0];
    randtype s1 = s[1]; //not constant, must be modificable var

    const randtype sum = s0 + s1;

    s1 ^= s0; //that's why s1 is not a const value.

    //rotation
    s[0] = rotl(s0,a) ^ s1 ^ (s1 << b);
    s[1] = rotl(s1,c);

    return sum;

}

double xorshiro::rand(double min, double max){

    return (double)(max-min) * ((double) xorshiro128p())
            / (double)(ULLONG_MAX) + min;
}
