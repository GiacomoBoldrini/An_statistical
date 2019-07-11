//c++ compute_and.cpp `root-config --libs --glibs --cflags` -o andamento

#include <iostream>
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TRandomGen.h"
#include "TLine.h"
#include "TF1.h"
#include "TLegend.h"
#include "xorshiro.h"

int main(){

    int samples = 100;
    double real_int = exp(1)-1;
    std::vector<double> deviazione;
    std::vector<double> error;
    std::vector<double> sam;
    
    xorshiro gen;
    
    for(int i = 0; i < 50; i++){

        double S_n = 0;
        double Error = 0;
        for(int j = 0; j < samples; j++){
            double x = gen.rand(0,1);
            S_n += exp(x);
            Error += pow(exp(x), 2);
        }

        S_n/=samples;

        double err_sn = sqrt(Error/samples - pow(S_n, 2))/sqrt(samples);

        double dev = S_n-real_int;
        deviazione.push_back(dev);
        std::cout<< "Integral: " << S_n << std::endl;
        std::cout << "Error: " << err_sn << std::endl;
        error.push_back(err_sn);
        sam.push_back(samples);

        samples += 50;
    }

    std::cout << sam.size() << " "<< deviazione.size() << " " << error.size() << std::endl;

    TGraphErrors* g = new TGraphErrors(sam.size(), &sam[0], &deviazione[0], 0, &error[0]);
    g->SetMarkerStyle(21);
    g->SetMarkerColor(4);
    g->GetYaxis()->SetTitleOffset(1.5);
    g->SetTitle("MC Estimator trend in function of the number of samples");
    g->GetXaxis()->SetTitle("Samples [N]");
    g->GetYaxis()->SetTitle("S_{n}-I");

    TGraph* g_error = new TGraph(sam.size(), &sam[0], &error[0]);
    g_error->SetMarkerStyle(21);
    g_error->SetMarkerColor(4);
    g_error->GetXaxis()->SetTitle("Samples [N]");
    g_error->GetYaxis()->SetTitle("Error");
    g_error->SetTitle("Error trend in function of the number of samples");

    TF1* sq = new TF1("#frac{1}{#sqrt{N}}", "[0]/sqrt(x[0])", 0, 2600);
    sq->SetLineColor(kRed);
    sq->SetLineStyle(2);

    TLine* real_val = new TLine(0, 0, 2600, 0);
    real_val->SetLineColor(kRed);
    real_val->SetLineStyle(2);

    TLegend* legend = new TLegend(.9,.9,0.6,0.7);
    legend->AddEntry(g_error, "Errors");
    legend->SetTextSize(.04);
    legend->AddEntry(sq, "#frac{1}{#sqrt{N}}");
    legend->SetTextSize(.04);

    g_error->Fit(sq);

    TLegend* legend1 = new TLegend(.9,.3,0.5,0.1);
    legend1->AddEntry(g, "MC Estimators");
    legend1->SetTextSize(.04);
    legend1->AddEntry(real_val, "Unbiased estimator");
    legend1->SetTextSize(.04);
    


    TCanvas*c = new TCanvas("c", "c", 1000,1000,1000,800);
    c->SetLogy();
    g_error->Draw("AP");
    //real_val->Draw();
    sq->Draw("same");
    legend->Draw();
    
    c->Update();
    c->Draw();
    c->SaveAs("./Error_trendtr.png");

    TCanvas*c1 = new TCanvas("c1", "c1", 1000,1000,1000,800);
    g->Draw("AP");
    real_val->Draw();
    legend1->Draw();
    c1->Draw();
    c1->SaveAs("./Estimator_trendtr.png");





    return 0;
}
