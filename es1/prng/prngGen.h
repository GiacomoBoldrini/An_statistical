#ifndef PRGGEN_
#define PRGGEN_
#include "prng.h"
//#include "xorshiro.h"


class prngGen{

    private:
    prng uniform; //oggetto della classe prng
    double pi;

    public:
    prngGen();
    ~prngGen();

    double BW(double x, double M0 = 0., double gamma =1.);

    double BW_inverse(double M0 = 0., double gamma = 1.);
    
    double BW_hm(double M0 = 0., double gamma = 1.,
                 double xmin = -10., double xmax = 10.,
                 double ymin = 0, double ymax = 0.1);

    double landau(double x, double mu = 0., double sigma = 1.);
    double landau_tc(double mu =0, double sigma = 1,
                          double xmin = 0, double xmax = 1,
                          double ymin = 0, double ymax = 0.5);

};

#endif
