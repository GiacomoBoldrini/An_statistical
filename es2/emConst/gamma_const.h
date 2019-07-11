#ifndef GAMMA_CONST_
#define GAMMA_CONST_

//#include "sum_alg.h"

class gamma_const{

    public:
        gamma_const();
        ~gamma_const();

        float gamma_constFI(int n = 100);
        double gamma_constDI (int n = 100);
        long double gamma_constLDI(int n = 100);
        float gamma_constFLLI(int n);
        double gamma_constDLLI(int n);
        long double gamma_constLDLLI(int n);
        float Kgamma_constFI(int n);
        long double Kgamma_constLDI(int n);
        float PWgamma_constFI(int n);
        long double PWgamma_constLDI(int n);
        
};




#endif 