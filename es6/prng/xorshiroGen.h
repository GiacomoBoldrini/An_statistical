#ifndef XORSHIROBWGEN_
#define XORSHIROBWGEN_
#include "prng.h"
#include "xorshiro.h"
#include <cmath>



class xorshiroGen{

    private:
    xorshiro uniform; //oggetto della classe xorshiro
    double pi;
    double* data ;
    double* smeared ;

    public:
    xorshiroGen();
    ~xorshiroGen();

    double BW(double x, double M0 = 0., double gamma =1.);

    double BW_inverse(double M0 = 0., double gamma = 1.);
    
    double BW_hm(double M0 = 0., double gamma = 1.,
                 double xmin = -10., double xmax = 10.,
                 double ymin = 0, double ymax = 0.1);

    double landau(double x, double mu = 0., double sigma = 1.);
    double landau_tc(double mu =0, double sigma = 1,
                          double xmin = 0, double xmax = 1,
                          double ymin = 0, double ymax = 0.5); 

    double gaus(double x, double mu=0, double sigma = 1);
    double gaus_tc(double mu =0, double sigma = 1, double xmin = -1, double xmax =1,
                    double ymin=1, double ymax=1);
    
    double decay(double x, double y, double alpha, double beta, double gamma );

    double gauss_bm(double mu, double sigma);


};

#endif
