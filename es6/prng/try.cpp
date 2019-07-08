#include "xorshiroGen.h"
#include <iostream>
#include <cmath>
#include "TH1D.h"
#include "TCanvas.h"
#include "simulation.h"


int main(int argc, char** argv){

    int N_ev = 50000;
    double diameter = 0.75e-6;
    double lambda = 400e-9; 
    double sm_fact = 2;

    simulation exper(N_ev,sm_fact,lambda,diameter);
    
    TCanvas* c = new TCanvas("c","c", 800,800);
    TH1D* h = new TH1D("h","h",100,0,M_PI);


    double* clean = exper.diffraction_tc(0,M_PI,0,20);
    double* smeared = exper.diffraction_smeared(0,M_PI,100);

    for(int i = 0; i<N_ev; i++){

        std::cout << clean[i] << std::endl;
        h->Fill(smeared[i]);
    }

    h->Draw("");
    c->Print("histo.pdf","pdf");




return 0;
}