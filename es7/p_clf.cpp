#include "percettrone.h"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "prng.h"
#include "xorshiro.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TLegend.h"
#include "TEllipse.h"
#include "TF1.h"

int main(){

    //----------------------------------
    //-------Retrieve data--------------
    //----------------------------------
    double x,y;

    TFile* dataset = new TFile("./dataset.root");
    TTree* signal = (TTree*) dataset->Get("Signal");
    signal->SetBranchAddress("x", &x);
    signal->SetBranchAddress("y", &y);
    TTree* background = (TTree*) dataset->Get("Background");
    background->SetBranchAddress("x", &x);
    background->SetBranchAddress("y", &y);

    int N_sig = signal->GetEntries();
    int N_bkg = background->GetEntries();


    //----------------------------------
    //-------Drawing result-------------
    //----------------------------------

    double xmin = -6;
    double xmax = 10;
    double ymin = -4;
    double ymax = 10;

    TCanvas* c1 = new TCanvas("c1","c1",1000,1000, 1000, 800);
    TH2F* gauss_sig = new TH2F("gauss_sig","gauss_sig", 100, xmin, xmax,100,ymin,ymax);
    gauss_sig->SetFillColor(kRed);
    TH2F* gauss_bkg = new TH2F("gauss_bkg","gauss_bkg",100,xmin,xmax,100,ymin,ymax);
    gauss_bkg->SetFillColor(kBlue);
    TLegend* legend = new TLegend(.1, .8, .4, .9);
    legend->AddEntry(gauss_sig, "Signal", "f");
    legend->AddEntry(gauss_bkg, "Background", "f");


    //----------------------------------
    //-----------Perceptron-------------
    //----------------------------------

    p_clf perceptron(2, N_sig+N_bkg);
    perceptron.set_weights();
    std::vector<double> w = perceptron.get_w();

    std::cout << w[0] << " " << w[1] << " " << w[2] << std::endl;

    double eta = 0.01;

    for(int i = 0; i < N_sig; i++){
        signal->GetEntry(i);
        perceptron.train(x, y, 0, eta);
        gauss_sig->Fill(x,y);

        background->GetEntry(i);
        perceptron.train(x, y, 1, eta);
        gauss_bkg->Fill(x,y);
        
    }

    w = perceptron.get_w();
    std::cout << w[0] << " " << w[1] << " " << w[2] << std::endl;

    //----------------------------------
    //-----------Accuracy---------------
    //----------------------------------

    int t_p = 0;
    int t_n = 0;
    int f_p = 0;
    int f_n = 0;

    for(int i = 0; i < N_sig; i++){

        signal->GetEntry(i);
        int p_sig = perceptron.predict(x, y);
        if (p_sig == 0){
            t_p += 1;
        }
        else{
            f_n += 1;
        }

        background->GetEntry(i);
        int p_bkg = perceptron.predict(x, y);
        if (p_bkg == 1){
            t_n += 1;
        }
        else{
            f_p += 1;
        }
        
    }

    std::cout << "TP: " << t_p << std::endl;
    std::cout << "FN: " << f_n << std::endl;
    std::cout << "TN: " << t_n << std::endl;
    std::cout << "FP: " << f_p << std::endl;
    std::cout<< "Accuracy: " << (double) (t_p+t_n)/(t_p+t_n+f_n+f_p) << std::endl;
    std::cout << "Error rate: " << (double) 1.- (t_p+t_n)/(t_p+t_n+f_n+f_p) << std::endl;


    //----------------------------------
    //----------Drawing DB--------------
    //----------------------------------

    line_bound line = perceptron.get_boundary();
    
    TF1* bound = new TF1("Decision boundary", "[0]*x + [1]", xmin, xmax);
    bound->SetParameters(line.m, line.q);
    bound->SetLineStyle(1);
    bound->SetLineWidth(2);
    bound->SetLineColor(kBlack);

    legend->AddEntry(bound, "Decision boundary", "l");

    
    gauss_sig->SetMarkerColor(kRed);
    gauss_bkg->SetMarkerColor(kBlue);
    gauss_sig->GetXaxis()->SetTitle("X [a.u.]");
    gauss_sig->GetYaxis()->SetTitle("Y [a.u.]");
    gauss_sig->SetTitle("Perceptron Decision Boundary");
    gauss_sig->SetStats(0000);
    gauss_sig->Draw();
    gauss_bkg->Draw("SAME");
    bound->Draw("same");
    legend->Draw();
    c1->Print("./analysis/perceptron.pdf","pdf");


    return 0;
}