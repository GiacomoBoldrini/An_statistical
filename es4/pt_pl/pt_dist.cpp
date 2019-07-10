//c++ pt_dist.cpp `root-config --libs --glibs --cflags` -o pt_dist

#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include <iostream>
#include <cmath>

double ptdist1(double *x, double *p){

    return (1/(5*M_PI))*log((sqrt(10*10-x[0]*x[0])+10)/abs(x[0]));
}

double ptdist2(double *x, double *p){

    return (1./10.)*log(10./x[0]);
}

int main(){

    TF1* pt1 = new TF1("pt_dist", ptdist1, 0, 10);
    TF1* pt2 = new TF1("pt_dist2", ptdist2, 0.01, 10.);
    pt1->SetLineColor(kBlue);
    pt1->SetLineWidth(2);
    pt2->SetLineColor(kRed);
    pt2->SetLineWidth(2);
    pt2->SetTitle("P_{T} PDF sampling difference; p_{T}; f_{P_{T}}(p_{T})");
    TCanvas* c = new TCanvas("c", "c", 1000,1000,1000,800);
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(pt1, "#theta ~ U[0, 2#pi]");
    legend->AddEntry(pt2, "sin(#theta) ~ U[-1,1]");
    pt2->Draw();
    pt1->Draw("same");
    
    legend->Draw();
    c->Draw();
    c->SaveAs("./pt_dist_2.png");

    pt1->SetTitle("P_{T} probability distribution function; p_{T}; f_{P_{T}}(p_{T})");
    pt1->Draw();
    c->Draw();
    c->SaveAs("./pt_dist.png");



    return 0;
}