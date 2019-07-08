#include <iostream>
#include "TF1.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <cmath>
#include "TLegend.h"

int main(){

    TF1* prima = new TF1("importance", "pow(x+1, 1.5)", 0, 1);
    prima->SetLineColor(kRed);
    prima->SetLineWidth(2);
    prima->SetTitle("Comparison between integrand and importance distribution");

    TF1* seconda = new TF1("Integrand", "exp(x)", 0, 1);
    seconda->SetLineColor(kBlue);
    seconda->SetLineWidth(2);

    TLegend* legend = new TLegend(.1,.9,0.4,0.7);
    legend->AddEntry(prima, "(x+1)^{1.5}");
    legend->SetTextSize(.04);
    legend->AddEntry(seconda, "e^{x}");
    legend->SetTextSize(.04);

    TCanvas* c = new TCanvas("c", "c", 1000,1000,1000,800);
    prima->Draw();
    prima->GetHistogram()->GetXaxis()->SetTitle("x");
    prima->GetHistogram()->GetYaxis()->SetTitle("f(x)");
    seconda->Draw("same");
    legend->Draw();
    c->Draw();
    c->SaveAs("./Importance2.png");
    return 0;
}