//c++ prng.cpp xorshiro.cpp xorshiroGen.cpp unfold.cpp `root-config --glibs --cflags` -L/home/wahid/Programs/RooUnfold/include -lRooUnfold -lUnfold -o unfold


#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>
#include "xorshiroGen.h"
#include <algorithm>
#include <vector>
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"


int main(int argc, char** argv){
    
    TFile* load = new TFile("/home/wahid/Wahid/University/analisi_statistica/esercizi/es6/data/data.root");

    TH1D* clean = (TH1D*)load->Get("clean");
    TH1D* smeared = (TH1D*)load->Get("smeared");

    int events = clean->GetEntries();

    RooUnfoldResponse response(smeared,clean);

    for(int i = 0; i<events; i++){

        std::cout << clean->GetBinContent(i) << std::endl;


    }
    



    return 0;
}