#include <climits>
#include "prng.h"
#include "prngGen.h"
#include "xorshiro.h"
#include "xorshiroGen.h"
#include <cmath>
#include "TMath.h"

xorshiroGen::xorshiroGen(){

    pi = atan(1)*4;

}

xorshiroGen::~xorshiroGen(){};

double xorshiroGen::BW_inverse(double M0, double gamma){

    return M0 + (gamma/2)*tan(pi*(uniform.rand() - 0.5));
}

double xorshiroGen::BW(double x, double M0, double gamma){

    return (gamma/2) /(pi*((x - M0)*(x - M0) + (gamma/2)*(gamma/2)));
}

double xorshiroGen::landau(double x, double mu, double sigma){


    return TMath::Landau(x,mu,sigma);
}

double xorshiroGen::BW_hm(double M0, double gamma,
                    double xmin , double xmax,
                    double ymin, double ymax){
                    
    double x_rand = uniform.rand(xmin,xmax);
    double y_rand = uniform.rand(ymin,ymax);
    double y = BW(x_rand, M0, gamma);
    
    while(y_rand > y){

        x_rand = uniform.rand(xmin,xmax);
        y_rand = uniform.rand(ymin,ymax);
        y = BW(x_rand,M0,gamma);

    }

    return x_rand;
}

double xorshiroGen::landau_tc(double mu, double sigma,
                          double xmin, double xmax,
                          double ymin, double ymax){

    double x_rand = uniform.rand(xmin,xmax);
    double y_rand = uniform.rand(ymin,ymax);
    double y = landau(x_rand, mu, sigma);

    while(y_rand > y){

        x_rand = uniform.rand(xmin,xmax);
        y_rand = uniform.rand(ymin,ymax);
        y = landau(x_rand,mu,sigma);
    }

    return x_rand;

}
    
    //double y_rand = (gamma/2) /((uniform.rand(xmin,xmax) - M0)*(uniform.rand(xmin,xmax)-M0) + (gamma/2)*(gamma/2));
