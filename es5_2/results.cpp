#include <iostream>
#include <cmath>
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"

int main(){

    std::vector<double> estim = {1.71835, 1.71829913, 1.71828458, 1.71826148};
    std::vector<double> estim_var = {0.00005, 0.00001013, 0.00000418, 0.00000880};

    std::vector<double> var = {0.24202856, 0.00212058, 0.000531207, 0.00391324};
    std::vector<double> err_var = {0.002296, 0.00002175, 0.00007748, 0.0005875};
    std::vector<double> number = {1,2,3,4};

    TLine* line = new TLine(1, exp(1)-1, 4, exp(1)-1);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);

    TGraphErrors* g = new TGraphErrors(4, &number[0], &estim[0], 0, &estim_var[0]);
    g->SetMarkerStyle(21);
    g->SetMarkerColor(4);
    g->SetTitle("MC estimators for different variance techniques");
    g->GetXaxis()->SetTitle("Variance reduction technique");
    g->GetYaxis()->SetTitle("Estimator Value");

    TGraphErrors* g_var = new TGraphErrors(4, &number[0], &var[0], 0, &err_var[0]);
    g_var->SetMarkerStyle(21);
    g_var->SetMarkerColor(4);
    g_var->SetTitle("MC estimator variance for different variance techniques");
    g_var->GetXaxis()->SetTitle("Variance reduction technique");
    g_var->GetYaxis()->SetTitle("Variance Value");

    TLegend* legend = new TLegend(.6, .5, .9, .9);
    legend->AddEntry(line, "#int_{0}^{1} e^{x}dx", "l");
    legend->AddEntry(g, "#splitline{1: Crude Mc}{#splitline{2: Stratified}{#splitline{3: Importance (x+1)^{1.5}}{4: Antithetic}}}");

    TCanvas* c = new TCanvas("c", "c", 1000,1000,1000,800);
    
    g->Draw("AP");
    line->Draw("same");
    legend->Draw();
    c->Draw();
    c->SaveAs("./Graphs/final_est_comp.png");

    c->SetLogy();
    g_var->Draw("AP");
    c->Draw();
    c->SaveAs("./Graphs/final_var_comp.png");

    return 0;
}