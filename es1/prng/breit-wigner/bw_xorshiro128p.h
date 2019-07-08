#ifndef XORSHIROBWGEN_
#define XORSHIROBWGEN_
#include "../prng.h"
#include "../xorshiro/xorshiro.h"


class xorshiroGen{

    private:
    xorshiro uniform; //oggetto della classe prng
    double pi;

    public:
    xorshiroGen();
    ~xorshiroGen();\

    double BW(double x, double M0 = 0., double gamma =1.);

    double BW_inverse(double M0 = 0., double gamma = 1.);
    
    double BW_hm(double M0 = 0., double gamma = 1.,
                 double xmin = -10., double xmax = 10.,
                 double ymin = 0, double ymax = 0.1);

};

#endif
