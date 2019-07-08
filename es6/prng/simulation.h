#ifndef SIMULATION_
#define SIMULATION_

#include "prng.h"
#include "xorshiro.h"
#include "xorshiroGen.h"
#include <cmath>
#include <iostream>

class simulation{
    public: 
        simulation(int N, double sm_fact_n,double lambda_n,double diameter_n);
        ~simulation();
        double diffraction_fun(double theta);
        double* diffraction_tc(double xmin,double xmax,double ymin,double ymax);

        double* diffraction_smeared(double xmin, double xmax, int bins);

    private:
        xorshiro uniform;
        xorshiroGen dist;
        int N_ev;
        double sm_fact;
        double lambda;
        double diameter;
        double* clean;
        double* smeared;
};



#endif //SIMULATION_