//c++ -o prova prova.cpp -L/usr/local/lib -lgsl -lgslcblas -lm
//c++ prng.cpp xorshiro.cpp xorshiroGen.cpp  prova.cpp -L/usr/local/lib -lgsl -lgslcblas -lm `root-config --libs --glibs --cflags` -o prova
// c++ prng.cpp xorshiro.cpp xorshiroGen.cpp  simulation.cpp main.cpp -L/usr/local/lib -lgsl -lgslcblas -lm `root-config --libs --glibs --cflags` -o prova
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
#include "TStyle.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "simulation.h"
#include <string>
#include <sstream>
#include "TPad.h"

TH1D* createRatio(TH1D* h1, TH1D* h2, const char* label="ratio"){
    TH1D* h3 = (TH1D*) h1->Clone("h3");
    h3->SetMarkerStyle(21);
    h3->SetMarkerColor(kBlack);
    h3->SetLineColor(kBlack);
    h3->SetTitle("");
    // Set up plot for markers and errors
    h3->Sumw2();
    h3->SetStats(0);
    h3->Divide(h2);
    //Fix automatically the ratio range
    h3->SetMaximum(h3->GetMaximum() + 0.05);
    h3->SetMinimum(h3->GetMinimum() - 0.05);
    
    // Adjust y-axis settings
    h3->GetYaxis()->SetTitle(label);
    return h3;
}

void line(){

    std::cout << "######################" << std::endl;
}

int main(int argc, char** argv){
 
 //------------------------------------------------------
 //-----------Generating distributions-------------------
 //------------------------------------------------------

//Simulation parameters
    int N = 50000;
    double diameter = 0.75e-6;
    double lambda = 400e-9;
    double k = (2*M_PI)/lambda;
    double c = 0.3;

    //Histo parameter
    double xmax = M_PI/2;
    double xmin = -M_PI/2;
    double ymax = 20;
    double ymin = 0;
    double bins = 80;

    double std = c * (xmax - xmin) / bins;

    if(argc > 1){

        N = std::stoi(argv[1]);
    }

    simulation experiment(N,std,lambda,diameter);

    double* clean = experiment.diffraction_tc(xmin,xmax,ymin,ymax);

    TH1D* h = new TH1D("Clean I(#theta)","Clean I(#theta)",bins,-M_PI/2,M_PI/2);
    h->SetTitle("I(#theta) Airy pattern");
    h->GetXaxis()->SetTitle("#theta");
    h->GetYaxis()->SetTitle("Counts/I_{0}");
    h->SetLineWidth(1);
    h->SetLineColor(kBlack);

    for(int i = 0; i < N; i++){

        h->Fill(clean[i]);
          
    }

    double m = h->GetMaximum();
    h->Scale(1./m);

    double* smeared = experiment.diffraction_smeared(); 

    
    TH1D* h_sm = new TH1D("histo_sm","histo_sm",bins,-M_PI/2,M_PI/2);
    h_sm->SetLineColor(kRed);
    for(int i = 0; i < N; i++){
        h_sm->Fill(smeared[i]);   
    }
    h_sm->Scale(1./m);

    experiment.set_smear(0.5*(xmax - xmin)/bins);
    double* smeared2 = experiment.diffraction_smeared(); 

    TH1D* h_sm2 = new TH1D("histo_sm2","histo_sm2",bins,-M_PI/2,M_PI/2);
    h_sm2->SetLineColor(kBlue);
    for(int i = 0; i < N; i++){
        h_sm2->Fill(smeared2[i]);   
    }
    h_sm2->Scale(1./m);

    experiment.set_smear(0.7*(xmax - xmin) / bins);
    double* smeared3 = experiment.diffraction_smeared(); 

    TH1D* h_sm3 = new TH1D("histo_sm3","histo_sm3",bins,-M_PI/2,M_PI/2);
    h_sm3->SetLineColor(kOrange);

    for(int i = 0; i < N; i++){
        h_sm3->Fill(smeared3[i]);   
    }

     h_sm3->Scale(1./m);

    experiment.set_smear(0.95*(xmax - xmin) / bins);
    double* smeared4 = experiment.diffraction_smeared(); 

    TH1D* h_sm4 = new TH1D("histo_sm4","histo_sm4",bins,-M_PI/2,M_PI/2);
    h_sm4->SetLineColor(kGreen);
    for(int i = 0; i < N; i++){
        h_sm4->Fill(smeared4[i]);   
    }

    h_sm4->Scale(1./m);
    

    //------------------------------------------------------
    //---------------------UNFOLDING------------------------
    //------------------------------------------------------

    RooUnfoldResponse response(bins,xmin,xmax); //Unfolding

    for(int i = 0; i < N; i++){
        response.Fill(smeared4[i],clean[i]);
    }

    //last par, iterations
    std::cout<< "Bayes" << std::endl;
    RooUnfoldBayes unfold(&response, h_sm4, 4);

    std::cout<< "BinbyBin" << std::endl;
    RooUnfoldBinByBin unfold_bbb(&response,h_sm4);

    std::cout<< "SVD" << std::endl;
    RooUnfoldSvd unfold_svd(&response,h_sm4);
    
    std::cout<< "------->UNFOLDING<-------" << std::endl;
    std::cout<< "------->Bayes<-------" << std::endl;
    TH1D* h_unfold =(TH1D*)unfold.Hreco();
    h_unfold->SetLineWidth(2);
    h_unfold->SetMarkerSize(2);
    std::cout<< "------->BBB<-------" << std::endl;
    TH1D* h_unfold_bbb = (TH1D*)unfold_bbb.Hreco();
    h_unfold_bbb->SetLineWidth(2);
    h_unfold_bbb->SetMarkerSize(2);
    std::cout<< "------->SVD<-------" << std::endl;
    //TH1D* h_unfold_svd = (TH1D*)unfold_svd.Hreco();
    //h_unfold_svd->SetLineWidth(2);
    //h_unfold_svd->SetMarkerSize(7);
    


    //------------------------------------------------------
    //------------------PLOTTING SMEARING-------------------
    //------------------------------------------------------

   TCanvas* c1 = new TCanvas("c1","c1",1000, 1000,1000, 800);
   //gStyle->SetOptStat(0000);
   TLegend* leg = new TLegend(0.1,0.7,0.4,0.9);

   leg->AddEntry(h, "Clean-data", "ple");
   leg->AddEntry(h_sm, "c = 0.3");
   leg->AddEntry(h_sm2, "c = 0.5");
   leg->AddEntry(h_sm3, "c = 0.7");
   leg->AddEntry(h_sm4, "c = 0.95");
   
   //Saving the results
   h->Draw();
   h_sm->Draw("hist same");
   h_sm2->Draw("hist same");
   h_sm3->Draw("hist same");
   h_sm4->Draw("hist same");
   //h->Write("clean");

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
   c1->SaveAs("./smeared.pdf");
   //data->Close();
    
    //------------------------------------------------------
    //------------------PLOTTING UNFOLDING------------------
    //------------------------------------------------------

    //TCanvas* c2 = new TCanvas ("c2", "c2", 1000,1000,1000,800);
    TLegend* leg1 = new TLegend(0.13,0.7,0.4,0.9);
    h_sm4->SetLineColor(kRed);
    h_sm4->SetLineWidth(2);
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    h_unfold->SetLineColor(kGreen);

    leg1->AddEntry(h, "True-data");
    leg1->AddEntry(h_sm4, "Smeared-data c = 0.95");
    leg1->AddEntry(h_unfold, "Bayes unfold", "ple");
    //leg1->AddEntry(h_unfold_bbb, "BinByBin unfold", "ple");
    //leg1->AddEntry(h_unfold_svd, "SVD unfold");

    TCanvas* c_bayes = new TCanvas("c_bayes", "c_bayes",1000, 1000,1000, 1000);
    // Upper histogram plot is pad1
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.28, 1, 1);
    pad1->SetBottomMargin(0);  // joins upper and lower plot
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.13);
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->Draw();
    // Lower ratio plot is pad2
    c_bayes->cd();  //returns to main canvas before defining pad2
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.28);
    pad2->SetTopMargin(0); // joins upper and lower plot
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.13);
    pad2->SetBottomMargin(0.13);
    pad2->SetGridx();
    pad2->SetGridy();
    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->Draw();

    c_bayes->Update();
    pad1->cd();

    gStyle->SetOptStat(0000);
    h->Draw("hist");
    h->SetTitle("Unfolding result using Bayesian approach");
    h->GetYaxis()->SetTitleSize(.04);
    h->GetXaxis()->SetTitleSize(.04);
    h_sm4->Draw("hist same");
    h_unfold->Draw("SAME");
    leg1->Draw();

    TH1D* rateo = createRatio(h, h_unfold);

    pad2->cd();
    rateo->Draw("hist");
    c_bayes->cd();
    c_bayes->Update();
    c_bayes->Draw();
    c_bayes->SaveAs("./provaBayes.pdf");
/* 
    h->SetLineColor(kBlue);
    h->Draw("hist");

    h_sm4->SetLineColor(kRed);
    h_sm4->Draw("hist same");

    h_unfold->SetLineColor(kGreen);
    h_unfold->Draw("SAME");

    h_unfold_bbb->SetLineColor(kPink);
    //h_unfold_bbb->Draw("SAME");

    //h_unfold_svd->SetLineColor(kOrange);
    //h_unfold_svd->Draw("hist SAME");
    //h_sm->Write("smeared");

    leg1->Draw("SAME");
    c2->Draw();
    c2->SaveAs("./Unfolding.pdf");

*/
    return 0;
     

     


}



 
