//c++ -o prova prova.cpp -L/usr/local/lib -lgsl -lgslcblas -lm
//c++ prng.cpp xorshiro.cpp xorshiroGen.cpp  prova.cpp -L/usr/local/lib -lgsl -lgslcblas -lm `root-config --libs --glibs --cflags` -o prova
#include <cmath>
#include <iostream>
#include "xorshiroGen.h"
#include <algorithm>
#include <vector>
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "simulation.h"
#include <string>
#include <sstream>

void line(){

    std::cout << "######################" << std::endl;
}

int main(int argc, char** argv){
 
//Simulation parameters
    int N = 50000;
    double diameter = 0.75e-6;
    double lambda = 400e-9;
    double c = 1;

    //Histo parameter
    double xmax = 3;
    double xmin = 0;
    double ymax = 20;
    double ymin = 0;
    double bins = 100;

    double std = c * (xmax - xmin) / bins;

    if(argc > 1){

        N = std::stoi(argv[1]);
    }

    xorshiroGen rndxor;
    xorshiro uniform;

    simulation experiment(N,c,lambda,diameter);

    TH1D* h = new TH1D("histo","histo",bins,xmin,xmax);
    TH1D* h_sm = new TH1D("histo_sm","histo_sm",bins,xmin,xmax);
    TCanvas* c1 = new TCanvas("c1","c1",700,700);

    double* clean = experiment.diffraction_tc(xmin,xmax,ymin,ymax);   
    double* smeared = experiment.diffraction_smeared(xmin,xmax,bins); 
    //double* diffraction_smeared(long int N,int bins,double xmin, double xmax,double sm_fact);

    RooUnfoldResponse response(bins,xmin,xmax); //Unfolding

    for(int i = 0; i < N; i++){

        double x = clean[i]; 

        h->Fill(clean[i]);
        h_sm->Fill(smeared[i]);
        response.Fill(smeared[i],clean[i]);     
    }
    //last par, iterations
    RooUnfoldBayes unfold(&response, h, 4);
    RooUnfoldBinByBin unfold_bbb(&response,h);
    RooUnfoldSvd unfold_svd(&response,h,15);

    TH1D* h_unfold =(TH1D*)unfold.Hreco();
    TH1D* h_unfold_bbb = (TH1D*)unfold_bbb.Hreco();
    TH1D* h_unfold_svd = (TH1D*)unfold_svd.Hreco();
    
   TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);

   leg->AddEntry(h, "Clean-data");
   leg->AddEntry(h_sm, "Smeared-data");
   leg->AddEntry(h_unfold, "Bayes unfold");
   leg->AddEntry(h_unfold_bbb, "BinByBin unfold");
   leg->AddEntry(h_unfold_svd, "SVD unfold");
   //Saving the results
   h_sm->SetLineColor(kRed);
   h_sm->Draw("");
   h->Draw("SAME");
   //h->Write("clean");
   
   h_unfold->SetLineColor(kGreen);   
   h_unfold->Draw("SAME");

   h_unfold_bbb->SetLineColor(kPink);
   h_unfold_bbb->Draw("SAME");

   h_unfold_svd->SetLineColor(kOrange);
   h_unfold_svd->Draw("SAME");
   //h_sm->Write("smeared");
   leg->Draw("SAME");
   c1->Draw("");

    std::ostringstream string_temp;
    string_temp << "_" << std << "-" 
        << lambda << "-" << diameter;
    std::string name_file = string_temp.str();

   std::cout << "\n\n Results simulation" << std::endl;
   line();
   line();
   std::cout << "Smearing factor: " << std << "\n" << "Lambda: " << lambda << "\n" << "Diameter: " << diameter << std::endl; 
   line();
   line();
   c1->SaveAs(("/home/wahid/Wahid/University/analisi_statistica/esercizi/es6/results/unfold"+name_file+".pdf").c_str(),"pdf");
   //data->Close();
   return 0;
     

     


}



 
