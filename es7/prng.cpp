#include "prng.h"
#include <climits>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


//Constructor 

prng::prng(){
        a = 6364136223846793005U;
        c = 1442695040888963407U;
        x = 1.;}


prng::~prng(){}
randtype prng::randInt(){
    x = (a * x + c);
    return x;
}

double prng::rand(double min, double max){

    return (double)(max-min) * ((double) randInt()) /
                ((double) ULLONG_MAX ) + min;

 }


