#include "monte_carlo.h"
#include <ctime>
#include <iostream>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLegend.h"
#include <cmath>

struct mean_pt{
    std::vector<double> pt;
    std::vector<double> std_err;
};

mean_pt generate_mean(std::vector<double> pt, std::vector<double> pl){
    
    mean_pt result;


    int bins = 100;
    double pl_min = 0;
    double pl_max = 10;
    double bin_step = (pl_max-pl_min)/bins;

    for(int i =0; i< bins-1; i++){
        std::vector<double> entries;
        int counter = 0;
        double bin_sum = 0;

        //Cycling on all vector
        for(int j=0; j<pt.size(); j++){
            
            if( (pl[j] > pl_min) & (pl[j] < pl_min+bin_step)){
                bin_sum += pt[j];
                entries.push_back(pt[j]);
                counter += 1;
                
            }
            
        }
        
        //computing variance of bin
        double var = 0;
        double mean = bin_sum/counter;
        for(int k = 0; k < entries.size(); k++){
            var+=(entries[k]-mean)*(entries[k]-mean);
        }
        double std_var = sqrt(var/(counter-1)/sqrt(counter));
        
        //pushing back in the vector
        
        result.pt.push_back(mean);
        result.std_err.push_back(std_var);
        pl_min += bin_step;
        
        
    }

    return result;

}

int main(){

    int N = 100000;
    std::vector<double> pt;
    std::vector<double> pl;

    TH2F* h2 = new TH2F("h2", "h2", 100, 0, 10, 100, 0, 10);
    h2->GetXaxis()->SetTitle("P_L [GeV/c]");
    h2->GetYaxis()->SetTitle("P_T [GeV/c]");
    h2->GetYaxis()->SetTitleOffset(2);
    h2->GetXaxis()->SetTitleOffset(2);
    h2->GetZaxis()->SetTitle("Counts");
    h2->GetZaxis()->SetTitleOffset(1.5);
    h2->SetTitle("P_T versus P_L distribution for N=100000 particles");

    //Equal generation but with sine uniform
    for(int i = 0; i < N; i++){
        generate_mc sim(rand());
        generation result = sim.simulation();
        h2->Fill(result.pl, result.pt);

    }

    TH2F* h3 = new TH2F("h3", "h3", 100, 0, 10, 100, 0, 10);
    h3->GetXaxis()->SetTitle("P_L [GeV/c]");
    h3->GetYaxis()->SetTitle("P_T [GeV/c]");
    h3->GetYaxis()->SetTitleOffset(2);
    h3->GetXaxis()->SetTitleOffset(2);
    h3->GetZaxis()->SetTitle("Counts");
    h3->GetZaxis()->SetTitleOffset(1.5);
    h3->SetTitle("P_T versus P_L distribution for N=100000 particles generating uniform sine");


    for(int i = 0; i < N; i++){
        generate_mc sim(rand());
        generation result = sim.simulation_sine();
        pt.push_back(result.pl);
        pl.push_back(result.pt);
        h3->Fill(result.pl, result.pt);

    }

    
    mean_pt averaging_pt = generate_mean(pt, pl);
    std::vector<double> pl_bins_err;
    std::vector<double> pl_binned;

    for(int i = 1; i < 100; i++){
        double bin_step = 10./100;
        pl_bins_err.push_back(bin_step);
        pl_binned.push_back(bin_step*i);

    }
    
    
    TGraphErrors* g = new TGraphErrors( 99, &pl_binned[0], &averaging_pt.pt[0], &pl_bins_err[0], &averaging_pt.std_err[0]);
    g->SetTitle("<P_{T}> Graph in funcion of P_{L} for 100000 MC events");
    g->GetXaxis()->SetTitle("P_{L} [GeV/c]");
    g->GetYaxis()->SetTitle("<P_{T}> [GeV/c]");
    g->SetMarkerColor(4);
    g->SetMarkerSize(1);
    g->SetMarkerStyle(20);
    
    TProfile *pt_m = h2->ProfileX("p_m", 0, 100);
    pt_m->GetXaxis()->SetTitle("P_{L} [GeV/c]");
    pt_m->GetXaxis()->SetTitleOffset(1.2);
    pt_m->GetYaxis()->SetTitle("<P_{T}> [GeV/c]");
    pt_m->SetTitle("<P_{T}> Graph in funcion of P_{L} for 100000 mc events");

    //TF1* analytical = new TF1("analytical", "((1-x*x)-x*TMath::ATan(sqrt((1/(x*x))-1)))/(-TMath::Log(x))", 0, 10);
    TF1* analytical = new TF1("analytical", "((1-x*x)-x*atan(sqrt(abs(1/(x*x)-1))))/-log(x)", 0, 10);
    analytical->SetLineWidth(2);
    analytical->SetLineColor(kRed);
    //TGraph* g = new TGraph( 99, &pl_binned[0], &averaging_pt.pt[0]);
    TCanvas*c = new TCanvas("c", "c", 1000, 1000, 1000, 800);
    //h2->Draw("colz");
    pt_m->Draw();
    //analytical->Draw();
    c->Draw();
    c->SaveAs("./p.png");

    gStyle->SetOptStat(0);
    c->SetLeftMargin(1);
    h2->Draw("lego");
    c->SaveAs("./gen_dist.png");

    // Sine generation
    TProfile *pt_m_s = h3->ProfileX("p_m_s", 0, 100);
    pt_m_s->GetXaxis()->SetTitle("P_{L} [GeV/c]");
    pt_m_s->GetXaxis()->SetTitleOffset(1.2);
    pt_m_s->GetYaxis()->SetTitle("<P_{T}> [GeV/c]");
    pt_m_s->SetTitle("<P_{T}> Graph in funcion of P_{L} for 100000 mc events uniform sine");
    pt_m_s->SetLineColor(kRed);
    pt_m->SetLineColor(kBlue);

    TLegend* l = new TLegend(0.7, 0.75, 0.9, 0.9);
    l->AddEntry(pt_m_s, "sin(#theta) ~ U[-1,1]");
    l->AddEntry(pt_m, "#theta ~ U[0, 2#pi]");
    pt_m->Draw("hist E");
    pt_m_s->Draw("hist E same");

    l->Draw();
    c->SaveAs("./gen_dist_size.png");

    h3->Draw("lego");
    c->SaveAs("./dist_sine.png");

    return 0;
}
