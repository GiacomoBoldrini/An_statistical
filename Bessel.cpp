#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <algorithm>
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "RooUnfoldResponse.h"

//#include "RooUnfoldResponse.h"



double bes(double*x, double* par){
    
    double diameter = 0.75e-6;
    double lambda = 400e-9;
    double k = (2*M_PI)/lambda;
    
    double alpha = k*diameter*sin(x[0]);
    double J_x = 2*gsl_sf_bessel_J1(alpha);
    
    return 1*((J_x)/alpha)*((J_x)/alpha);
    
}

double smear(double x_true, double std){
    double x_smear = gRandom->Gaus(0,std);
    return x_true+x_smear;
}

int main(){
    
    std::vector<double> x_;
    std::vector<double> y_;
    
    /*
    double h = -M_PI/2;
    double step = M_PI/3000;
    
    for(int i =0; i < 3000; i++){
        
        y_.push_back(bes(h));
        x_.push_back(h);
        h+=step;
        
    }
    
    TCanvas* c = new TCanvas("c", "c", 1000, 1000, 1000, 800);
    TGraph* g = new TGraph(x_.size(), &x_[0], &y_[0]);
    TLine* line = new TLine(0, 0, 30, 0);
    line->SetLineStyle(2);
    
     */
    TF1* fun = new TF1("f1", &bes, -M_PI/2, M_PI/2, 0);
    TH1F* h = new TH1F("h", "h", 80, -M_PI/2, M_PI/2);
    h->SetLineColor(kBlue);
    TH1F* h_sm = new TH1F("h_s", "h_s", 80, -M_PI/2, M_PI/2);
    h_sm->SetLineColor(kRed);
    TH1F* h_sm1 = new TH1F("h_s1", "h_s1", 80, -M_PI/2, M_PI/2);
    h_sm1->SetLineColor(kOrange);
    TH1F* h_sm2 = new TH1F("h_s2", "h_s2", 80, -M_PI/2, M_PI/2);
    h_sm2->SetLineColor(kGreen);
    TH1F* h_sm3 = new TH1F("h_s3", "h_s3", 80, -M_PI/2, M_PI/2);
    h_sm3->SetLineColor(kPink);
    
    double max = 0;
    double min = 0;
    double bin_wid = 0.0389846;
    
    
    for(int i = 0; i < 10000; i++){
        double x = fun->GetRandom();
        if(x > max) max = x;
        if(x < min) min = x;
        h->Fill(x);
        h_sm->Fill(smear(x, 0.4*bin_wid));
        h_sm1->Fill(smear(x, 0.5*bin_wid));
        h_sm2->Fill(smear(x, 0.6*bin_wid));
        h_sm3->Fill(smear(x, 0.7*bin_wid));
    }
    
    
    std::cout << (max-min)/80 <<std::endl;
    TCanvas* c = new TCanvas("c", "c", 1000, 1000, 1000, 800);
    h->Draw("hist");
    h_sm->Draw("hist same");
    h_sm1->Draw("hist same");
    h_sm2->Draw("hist same");
    h_sm3->Draw("hist same");
    c->Draw();
    c->SaveAs("./prova.png");
//
//
//    g->SetMarkerStyle(20);
//    g->SetMarkerColor(4);
//    g->SetMarkerSize(0.5);
//    g->SetTitle("Theoretical I(#theta) Airy function; #theta; I(#theta)");
//
//    g->Draw("AP");
//    //line->Draw("same");
//    c->Draw();
//    c->SaveAs("./Tairy.pdf");
//
    RooUnfoldResponse* responde= new RooUnfoldResponse(80,-M_PI/2, M_PI/2);
    
    return 0;
    
}
