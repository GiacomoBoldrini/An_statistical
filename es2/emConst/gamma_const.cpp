#include "sum_alg.h"
#include "gamma_const.h"
#include <cmath>
#include <iostream>
#include <math.h>


gamma_const::gamma_const(){}
gamma_const::~gamma_const(){}

double gamma_const::gamma_constDI(int n){

    sum <double,int> sumDI;
    return  sumDI.norm_sum(n) - log(n);;
}

float gamma_const::gamma_constFI(int n){

    sum <float,int> sumDI;
    return  sumDI.norm_sum(n) - log(n);
}

long double gamma_const::gamma_constLDI(int n){

    sum <long double,int> sumDI;
    return  sumDI.norm_sum(n) - log(n);;
}

float gamma_const::gamma_constFLLI(int n){

    sum <float,long long int> sumDI;
    return  sumDI.norm_sum(n) - log(n);
}

double gamma_const::gamma_constDLLI(int n){

    sum <double,long long int> sumDI;
    return  sumDI.norm_sum(n) - log(n);
}

long double gamma_const::gamma_constLDLLI(int n){

    sum <long double,long long int> sumDI;
    return  sumDI.norm_sum(n) - log(n);
}

long double gamma_const::Kgamma_constLDI(int n){

    sum <long double,long long int> sumDI;
    return  sumDI.kahan_sum(n) - log(n);
}

float gamma_const::Kgamma_constFI(int n){

    sum <float,long long int> sumDI;
    return  sumDI.kahan_sum(n) - log(n);
}

float gamma_const::PWgamma_constFI(int n){

    sum <float,long long int> sumDI;

    return sumDI.pairwise_sum(1,n) - log(n);
}

long double gamma_const::PWgamma_constLDI(int n){

    sum <long double,long long int> sumDI;

    return sumDI.pairwise_sum(1,n) - log(n);
}



