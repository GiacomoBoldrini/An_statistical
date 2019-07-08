#include <climits>
#include "prng.h"
//#include "prngGen.h"
#include "xorshiro.h"
#include "xorshiroGen.h"
#include <cmath>
#include "TMath.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <algorithm>
#include "TRandom3.h"
#include <ctime>


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

double xorshiroGen::gaus(double x, double mu, double std){

    double d = x - mu;

    return (1/sqrt(2*M_PI*std)*exp(-d*d/(2*std*std)));
}

double xorshiroGen::gauss_bm(double mu, double sigma){
    // box-muller algorithm
    double u1 = uniform.rand(0,1);
    double u2 = uniform.rand(0,1);
    double z0 = sqrt( -2.0 * log(u1) ) * cos (2*M_PI*u2);
    // double z1: discarded
    return z0 * sigma + mu ;
}


double xorshiroGen::gaus_tc(double mu,double std,double xmin,double xmax,double ymin, double ymax){

    double x = 0;
    double y = 0;
    

    do
    {
        x = uniform.rand(xmin,xmax);
        y = uniform.rand(ymin,ymax);
    } while (y > gaus(x,mu,std));
    
    return x;


}

double xorshiroGen::decay(double theta, double phi, double alpha, double beta, double gamma ){

    return (3./(4.*M_PI))*(0.5*(1.-alpha) + (0.5)*(3.*alpha-1)*cos(theta)*cos(theta) - beta*sin(theta)*sin(theta)*cos(2.*phi)- sqrt(2.)*gamma*sin(2.*theta)*cos(phi));

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
