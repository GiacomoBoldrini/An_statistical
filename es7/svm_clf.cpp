#include "svm.h"
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
#include "TF2.h"
#include <array>

double svm_hyp(std::vector<double> w, double x, double y){
        return w[0] + sqrt(2)*w[1]*x + sqrt(2)*w[2]*y + w[3]*x*x + sqrt(2)*w[4]*x*y+w[5]*y*y;
    }


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
//-----------SVM-------------
//----------------------------------

    svm_clf svm(N_sig+N_bkg, 2, 10000000);
    svm.set_weights();
    std::vector<double> w = svm.get_weight();

    std::cout << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << " " << w[4] << " " << w[5] << std::endl;
    double eta = 0.01;

    double accuracy = 0;

    int counter = 0;
    while(accuracy != 1){

        if(counter == 10){
            svm.momentum();
            counter = 0;
            std::cout << "Momentum" <<std::endl;
        }
        for(int j = 0; j < 10; j++){
            for(int i = 0; i < N_sig; i++){
                signal->GetEntry(i);
                svm.train(x, y, 1, eta);
                if(j == 0){
                    gauss_sig->Fill(x,y);
                }

                background->GetEntry(i);
                svm.train(x, y, -1, eta);
                if(j == 0){
                    gauss_bkg->Fill(x,y);
                }
                
            }
        }

        w = svm.get_weight();

        std::cout << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << " " << w[4] << " " << w[5] << std::endl;

        double t_p = 0.;
        double t_n = 0.;
        double f_p = 0.;
        double f_n = 0.;

        for(int i = 0; i < N_sig; i++){

            signal->GetEntry(i);
            int p_sig = svm.predict(x, y);
            if (p_sig == 1){
                t_p += 1;
            }
            else{
                f_n += 1;
                std::cout << x << " " << y << std::endl;
            }

            background->GetEntry(i);
            int p_bkg = svm.predict(x, y);
            if (p_bkg == -1){
                t_n += 1;
            }
            else{
                f_p += 1;
                std::cout << x << " " << y << std::endl;
            }
            
        }
        accuracy = (t_p+t_n)/(t_p+t_n+f_n+f_p);
        counter ++;
        std::cout << accuracy << std::endl;

    }


//----------------------------------
//-----------Accuracy---------------
//----------------------------------
/* 
    double t_p = 0.;
    double t_n = 0.;
    double f_p = 0.;
    double f_n = 0.;

    for(int i = 0; i < N_sig; i++){

        signal->GetEntry(i);
        int p_sig = svm.predict(x, y);
        if (p_sig == 1){
            t_p += 1;
        }
        else{
            f_n += 1;
        }

        background->GetEntry(i);
        int p_bkg = svm.predict(x, y);
        if (p_bkg == -1){
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
    double accuracy = (t_p+t_n)/(t_p+t_n+f_n+f_p);
    double error = 1.-accuracy;
    std::cout<< "Accuracy: " << accuracy << std::endl;
    std::cout << "Error rate: " << error << std::endl;

    

*/
    w = svm.get_weight();
    TF2* db = new TF2("Decision Boundary SVM", "[0] + sqrt(2)*[1]*x + sqrt(2)*[2]*y + [3]*x*x + sqrt(2)*[4]*x*y+[5]*y*y", xmin, xmax, ymin, ymax);
    db->SetParameters(w[0], w[1], w[2], w[3], w[4], w[5]);
    std::vector<double> level {5.};

    db->SetContour(100);

    
    db->SetContour(4);
    gauss_sig->SetMarkerColor(kRed);
    gauss_bkg->SetMarkerColor(kBlue);
    gauss_sig->GetXaxis()->SetTitle("X [a.u.]");
    gauss_sig->GetYaxis()->SetTitle("Y [a.u.]");
    gauss_sig->SetTitle("Perceptron Decision Boundary");
    gauss_sig->SetStats(0000);

    TCanvas* c = new TCanvas("c", "c", 1000,1000,1000, 800);
    
    gauss_sig->Draw();
    gauss_bkg->Draw("SAME");
    db->Draw("same");
    legend->Draw();
    c->Draw();
    c->SaveAs("./prova.png");
    return 0;
}