//c++ results.cpp `root-config --libs --glibs --cflags` -o res

#include <iostream>
#include <cmath>
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"

int main(){

    std::vector<double> estim = {1.71820764, 1.71826621, 1.71829322, 1.71828255};
    std::vector<double> estim_var = {0.00011599, 0.00001030, 0.00000430, 0.00000888};

    std::vector<double> var = {0.24198802, 0.00211949, 0.00053135, 0.00391174};
    std::vector<double> err_var = {0.00002338, 0.00000022, 0.00000008, 0.00000060};
    std::vector<double> number = {1,2,3,4};

    TLine* line = new TLine(1, exp(1)-1, 4, exp(1)-1);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);

    TGraphErrors* g = new TGraphErrors(4, &number[0], &estim[0], 0, &estim_var[0]);
    g->SetMarkerStyle(21);
    g->SetMarkerColor(4);
    g->SetTitle("MC estimators for different variance techniques xorshiro128p algorithm");
    g->GetXaxis()->SetTitle("Variance reduction technique");
    g->GetYaxis()->SetTitle("Estimator Value");

    TGraphErrors* g_var = new TGraphErrors(4, &number[0], &var[0], 0, &err_var[0]);
    g_var->SetMarkerStyle(21);
    g_var->SetMarkerColor(4);
    g_var->SetTitle("MC estimator variance for different variance techniques xorshiro128p algorithm");
    g_var->GetXaxis()->SetTitle("Variance reduction technique");
    g_var->GetYaxis()->SetTitle("Variance Value");

    TLegend* legend = new TLegend(.5, .1, .9, .5);
    legend->AddEntry(line, "#int_{0}^{1} e^{x}dx", "l");
    legend->AddEntry(g, "#splitline{1: Crude Mc}{#splitline{2: Stratified}{#splitline{3: Importance (x+1)^{1.5}}{4: Antithetic}}}");

    TCanvas* c = new TCanvas("c", "c", 1000,1000,1000,800);
    
    g->Draw("AP");
    line->Draw("same");
    legend->Draw();
    c->Draw();
    c->SaveAs("./Graphs/final_est_comp_xor.png");

    c->SetLogy();
    g_var->Draw("AP");
    c->Draw();
    c->SaveAs("./Graphs/final_var_comp_xor.png");

    return 0;
}
