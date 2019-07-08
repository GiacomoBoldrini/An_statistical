#ifndef XOROSHIRO_
#define XOROSHIRO_

#include "prng.h"

typedef long long unsigned int randtype;

class xorshiro
{
private:

    prng rnd;
    int a;
    int b;
    int c;
    randtype s[2];
    


public:
    xorshiro();
    ~xorshiro();

    randtype rotl(const randtype, int k);
    randtype xorshiro128p();

    double rand(double min = 0., double max = 1.);   

    double landau(double x, double mu = 0., double sigma = 1.);
    double landau_tc(double mu =0, double sigma = 1,
                          double xmin = 0, double xmax = 1,
                          double ymin = 0, double ymax = 0.5); 


};

#endif //XOROSHIRO_
