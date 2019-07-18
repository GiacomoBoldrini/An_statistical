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
//#include "RooUnfoldIds.h"
#include "simulation.h"
#include <string>
#include <sstream>
#include "TPad.h"

TH1D* createRatio(TH1D* h1, TH1D* h2, const char* label="ratio"){
    TH1D* h3 = (TH1D*) h1->Clone("h3");
    h3->SetMarkerStyle(21);
    h3->SetMarkerSize(1);
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
    simulation experiment_dummy(N,0.1,lambda,diameter);

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

    TH1D* Intensity_measured = new TH1D("Intensity Measured", "Intensity Measured", bins, xmin, xmax);

    experiment.set_smear(0.95*(xmax - xmin)/bins);
    double* smeared_to_response = experiment.diffraction_smeared();
    double* smeared_to_unfold = experiment.diffraction_smeared();

    for(int i = 0; i < N; i++){

        Intensity_measured->Fill(smeared_to_unfold[i]);
        response.Fill(smeared_to_response[i],clean[i]);
    }

    Intensity_measured->Scale(1./m);

    //last par, iterations
    std::cout<< "Bayes" << std::endl;
    RooUnfoldBayes unfold(&response, Intensity_measured, 4);

    std::cout<< "BinbyBin" << std::endl;
    RooUnfoldBinByBin unfold_bbb(&response,Intensity_measured);

    std::cout<< "SVD" << std::endl;
    RooUnfoldSvd unfold_svd(&response,Intensity_measured, 20);

    //std::cout<< "IDS" << std::endl;
    //RooUnfoldIds unfold_ids(&response,Intensity_measured, 4);
    
    std::cout<< "------->UNFOLDING<-------" << std::endl;
    std::cout<< "------->Bayes<-------" << std::endl;
    TH1D* h_unfold =(TH1D*)unfold.Hreco();
    h_unfold->SetLineWidth(2);
    h_unfold->SetMarkerSize(2);
    std::cout<< "------->BBB<-------" << std::endl;
    TH1D* h_unfold_bbb = (TH1D*)unfold_bbb.Hreco();
    h_unfold_bbb->SetLineWidth(2);
    h_unfold_bbb->SetMarkerSize(2);
    //std::cout<< "------->SVD<-------" << std::endl;
    //TH1D* h_unfold_svd = (TH1D*)unfold_svd.Hreco();
    //h_unfold_svd->SetLineWidth(2);
    //h_unfold_svd->SetMarkerSize(7);
    //std::cout<< "------->IDS<-------" << std::endl;
    //TH1D* h_unfold_ids = (TH1D*)unfold_ids.Hreco();
    //h_unfold_ids->SetLineWidth(2);
    //h_unfold_ids->SetMarkerSize(7);
    


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
    Intensity_measured->SetLineColor(kRed);
    Intensity_measured->SetLineWidth(1);
    h->SetLineColor(kBlue);
    h->SetLineWidth(1);
    h_unfold->SetLineColor(kGreen);

    leg1->AddEntry(h, "True-data");
    leg1->AddEntry(Intensity_measured, "Smeared-data c = 0.95");
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

    gStyle->SetOptStat(0);
    h->Draw("hist");
    h->SetTitle("Unfolding result using Bayesian approach");
    h->GetYaxis()->SetTitleSize(.04);
    h->GetXaxis()->SetTitleSize(.04);
    Intensity_measured->Draw("hist same");
    h_unfold->Draw("SAME");
    leg1->Draw();

    TH1D* rateo = createRatio(h, h_unfold);

    pad2->cd();
    rateo->Draw("hist");
    TH1F* hratioerror =(TH1F*) rateo->DrawCopy("E2 same");
    hratioerror->SetFillStyle(3013);
    hratioerror->SetFillColor(13);
    hratioerror->SetMarkerStyle(1);
    c_bayes->cd();
    c_bayes->Update();
    c_bayes->Draw();
    c_bayes->SaveAs("./provaBayes.pdf");



    //---------------------------------
    //---------Unfolding BinbyBin------
    //---------------------------------

    TLegend* leg2 = new TLegend(0.13,0.7,0.4,0.9);
    Intensity_measured->SetLineColor(kRed);
    Intensity_measured->SetLineWidth(1);
    h->SetLineColor(kBlue);
    h->SetLineWidth(1);
    h_unfold_bbb->SetLineColor(kGreen);

    leg2->AddEntry(h, "True-data");
    leg2->AddEntry(Intensity_measured, "Smeared-data c = 0.95");
    leg2->AddEntry(h_unfold_bbb, "Bin-by_Bin unfold", "ple");
    //leg1->AddEntry(h_unfold_bbb, "BinByBin unfold", "ple");
    //leg1->AddEntry(h_unfold_svd, "SVD unfold");

    TCanvas* c_bbb = new TCanvas("c_bbb", "c_bbb",1000, 1000,1000, 1000);
    // Upper histogram plot is pad1
    //TPad* pad1 = new TPad("pad1", "pad1", 0, 0.28, 1, 1);
    pad1->SetBottomMargin(0);  // joins upper and lower plot
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.13);
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->Draw();
    // Lower ratio plot is pad2
    c_bbb->cd();  //returns to main canvas before defining pad2
    //TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.28);
    pad2->SetTopMargin(0); // joins upper and lower plot
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.13);
    pad2->SetBottomMargin(0.13);
    pad2->SetGridx();
    pad2->SetGridy();
    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->Draw();

    c_bbb->Update();
    pad1->cd();

    gStyle->SetOptStat(0);
    h->Draw("hist");
    h->SetTitle("Unfolding result using Bin-by-Bin approach");
    h->GetYaxis()->SetTitleSize(.04);
    h->GetXaxis()->SetTitleSize(.04);
    Intensity_measured->Draw("hist same");
    h_unfold_bbb->Draw("SAME");
    leg2->Draw();

    TH1D* rateo_b = createRatio(h, h_unfold_bbb);

    pad2->cd();
    rateo_b->Draw("hist");
    TH1F* hratioerror_b =(TH1F*) rateo_b->DrawCopy("E2 same");
    hratioerror_b->SetFillStyle(3013);
    hratioerror_b->SetFillColor(13);
    hratioerror_b->SetMarkerStyle(1);
    c_bbb->cd();
    c_bbb->Update();
    c_bbb->Draw();
    c_bbb->SaveAs("./provaBBB.pdf");

/* 
    //-----------------------------------
    //-------Plotting SVD----------------
    //-----------------------------------


    TLegend* leg3 = new TLegend(0.13,0.7,0.4,0.9);
    Intensity_measured->SetLineColor(kRed);
    Intensity_measured->SetLineWidth(1);
    h->SetLineColor(kBlue);
    h->SetLineWidth(1);
    h_unfold_svd->SetLineColor(kGreen);

    leg3->AddEntry(h, "True-data");
    leg3->AddEntry(Intensity_measured, "Smeared-data c = 0.95");
    leg3->AddEntry(h_unfold_svd, "SVD unfold", "ple");
    //leg1->AddEntry(h_unfold_bbb, "BinByBin unfold", "ple");
    //leg1->AddEntry(h_unfold_svd, "SVD unfold");

    TCanvas* c_svd = new TCanvas("c_svd", "c_svd",1000, 1000,1000, 1000);
    // Upper histogram plot is pad1
    //TPad* pad1 = new TPad("pad1", "pad1", 0, 0.28, 1, 1);
    pad1->SetBottomMargin(0);  // joins upper and lower plot
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.13);
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->Draw();
    // Lower ratio plot is pad2
    c_svd->cd();  //returns to main canvas before defining pad2
    //TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.28);
    pad2->SetTopMargin(0); // joins upper and lower plot
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.13);
    pad2->SetBottomMargin(0.13);
    pad2->SetGridx();
    pad2->SetGridy();
    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->Draw();

    c_svd->Update();
    pad1->cd();

    gStyle->SetOptStat(0);
    h->Draw("hist");
    h->SetTitle("Unfolding result using SVD approach");
    h->GetYaxis()->SetTitleSize(.04);
    h->GetXaxis()->SetTitleSize(.04);
    Intensity_measured->Draw("hist same");
    h_unfold_svd->Draw("SAME");
    leg3->Draw();

    TH1D* rateo_svd = createRatio(h, h_unfold_svd);

    pad2->cd();
    rateo_svd->Draw("hist");
    TH1F* hratioerror_svd =(TH1F*) rateo_svd->DrawCopy("E2 same");
    hratioerror_svd->SetFillStyle(3013);
    hratioerror_svd->SetFillColor(13);
    hratioerror_svd->SetMarkerStyle(1);
    c_svd->cd();
    c_svd->Update();
    c_svd->Draw();
    c_svd->SaveAs("./provaSVD.pdf");


 

    //-----------------------------------
    //-------Plotting IDS----------------
    //-----------------------------------


    TLegend* leg4 = new TLegend(0.13,0.7,0.4,0.9);
    Intensity_measured->SetLineColor(kRed);
    Intensity_measured->SetLineWidth(1);
    h->SetLineColor(kBlue);
    h->SetLineWidth(1);
    h_unfold_ids->SetLineColor(kGreen);

    leg4->AddEntry(h, "True-data");
    leg4->AddEntry(Intensity_measured, "Smeared-data c = 0.95");
    leg4->AddEntry(h_unfold_ids, "IDS unfold", "ple");
    //leg1->AddEntry(h_unfold_bbb, "BinByBin unfold", "ple");
    //leg1->AddEntry(h_unfold_svd, "SVD unfold");

    TCanvas* c_ids = new TCanvas("c_ids", "c_ids",1000, 1000,1000, 1000);
    // Upper histogram plot is pad1
    //TPad* pad1 = new TPad("pad1", "pad1", 0, 0.28, 1, 1);
    pad1->SetBottomMargin(0);  // joins upper and lower plot
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.13);
    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->Draw();
    // Lower ratio plot is pad2
    c_ids->cd();  //returns to main canvas before defining pad2
    //TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.28);
    pad2->SetTopMargin(0); // joins upper and lower plot
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.13);
    pad2->SetBottomMargin(0.13);
    pad2->SetGridx();
    pad2->SetGridy();
    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->Draw();

    c_ids->Update();
    pad1->cd();

    gStyle->SetOptStat(0);
    h->Draw("hist");
    h->SetTitle("Unfolding result using IDS approach");
    h->GetYaxis()->SetTitleSize(.04);
    h->GetXaxis()->SetTitleSize(.04);
    Intensity_measured->Draw("hist same");
    h_unfold_ids->Draw("SAME");
    leg4->Draw();

    TH1D* rateo_ids = createRatio(h, h_unfold_ids);

    pad2->cd();
    rateo_ids->Draw("hist");
    TH1F* hratioerror_ids =(TH1F*) rateo_ids->DrawCopy("E2 same");
    hratioerror_ids->SetFillStyle(3013);
    hratioerror_ids->SetFillColor(13);
    hratioerror_ids->SetMarkerStyle(1);
    c_ids->cd();
    c_ids->Update();
    c_ids->Draw();
    c_ids->SaveAs("./provaIDS.pdf");

*/

    //--------------------------------
    //-------PLOT COVARIANCE----------
    //--------------------------------

    //Bayes
    TCanvas* c_cov_bayes = new TCanvas("c_c_b", "c_c_b",1000, 1000,1000, 1000);
    unfold.Ereco().Draw("colz");
    c_cov_bayes->Draw();
    c_cov_bayes->SaveAs("./cov_bayes.pdf");

    unfold.ErecoV().Draw("colz");
    c_cov_bayes->Draw();
    c_cov_bayes->SaveAs("./cov_bayes_V.pdf");

    //Bin-to-Bin
    TCanvas* c_cov_bbb = new TCanvas("c_c_bbb", "c_c_bbb",1000, 1000,1000, 1000);
    unfold_bbb.Ereco().Draw("colz");
    c_cov_bbb->Draw();
    c_cov_bbb->SaveAs("./cov_bbb.pdf");

    unfold_bbb.ErecoV().Draw("colz");
    c_cov_bbb->Draw();
    c_cov_bbb->SaveAs("./cov_bbb_V.pdf");
/* 
    //SVD
    TCanvas* c_cov_svd = new TCanvas("c_c_svd", "c_c_svd",1000, 1000,1000, 1000);
    unfold_svd.Ereco().Draw("colz");
    c_cov_svd->Draw();
    c_cov_svd->SaveAs("./cov_svd.pdf");

    unfold_svd.ErecoV().Draw("colz");
    c_cov_svd->Draw();
    c_cov_svd->SaveAs("./cov_svd_V.pdf");


    //IDS
    TCanvas* c_cov_ids = new TCanvas("c_c_ids", "c_c_ids",1000, 1000,1000, 1000);
    unfold_ids.Ereco().Draw("colz");
    c_cov_ids->Draw();
    c_cov_ids->SaveAs("./cov_ids.pdf");

    unfold_ids.ErecoV().Draw("colz");
    c_cov_ids->Draw();
    c_cov_ids->SaveAs("./cov_ids_V.pdf");
*/

    //------------------------------
    //------Hypervariation of c-----
    //------------------------------


    RooUnfoldResponse response1(bins,xmin,xmax); //Unfolding
    RooUnfoldResponse response2(bins,xmin,xmax); //Unfolding
    RooUnfoldResponse response3(bins,xmin,xmax); //Unfolding

    TH1D* Intensity_measured1 = new TH1D("Intensity Measured1", "Intensity Measured1", bins, xmin, xmax);
    TH1D* Intensity_measured2 = new TH1D("Intensity Measured2", "Intensity Measured2", bins, xmin, xmax);
    TH1D* Intensity_measured3 = new TH1D("Intensity Measured3", "Intensity Measured3", bins, xmin, xmax);

    double* clean_dummy = experiment_dummy.diffraction_tc(xmin,xmax,ymin,ymax);
    double* smeared_conf = experiment_dummy.diffraction_smeared();

    //-------0.1------------
    experiment.set_smear(0.1);
    double* smeared_to_response1 = experiment.diffraction_smeared();
    //double* smeared_to_unfold1 = experiment.diffraction_smeared();

    for(int i = 0; i < N; i++){

        Intensity_measured1->Fill(smeared_conf[i]);
        response1.Fill(smeared_to_response1[i],clean[i]);
    }

    Intensity_measured1->Scale(1./m);

    std::cout<< "Bayes" << std::endl;
    RooUnfoldBayes unfold1(&response1, Intensity_measured1, 4);

    //-------0.3------------
    experiment.set_smear(0.3);
    experiment_dummy.set_smear(0.3);
    double* smeared_to_response2 = experiment.diffraction_smeared();
    //double* smeared_to_unfold2 = experiment.diffraction_smeared();
    double* smeared_to_unfold2 = experiment_dummy.diffraction_smeared();

    for(int i = 0; i < N; i++){

        Intensity_measured2->Fill(smeared_to_unfold2[i]);
        response2.Fill(smeared_to_response2[i],clean[i]);
    }

    Intensity_measured2->Scale(1./m);

    std::cout<< "Bayes" << std::endl;
    RooUnfoldBayes unfold2(&response2, Intensity_measured2, 4);


    //-------0.5------------
    experiment.set_smear(0.5);
    experiment_dummy.set_smear(0.5);
    double* smeared_to_response3 = experiment.diffraction_smeared();
    //double* smeared_to_unfold3 = experiment.diffraction_smeared();
    double* smeared_to_unfold3 = experiment_dummy.diffraction_smeared();

    for(int i = 0; i < N; i++){

        Intensity_measured3->Fill(smeared_to_unfold3[i]);
        response3.Fill(smeared_to_response3[i],clean[i]);
    }

    Intensity_measured3->Scale(1./m);

    std::cout<< "Bayes" << std::endl;
    RooUnfoldBayes unfold3(&response3, Intensity_measured3, 4);


    //---Unfolding------------

    std::cout << "HyperUnfolding " << std::endl;

    TH1D* h_unfold1 =(TH1D*)unfold1.Hreco();
    h_unfold1->SetLineWidth(2);
    h_unfold1->SetMarkerSize(2);

    TH1D* h_unfold2 =(TH1D*)unfold2.Hreco();
    h_unfold2->SetLineWidth(2);
    h_unfold2->SetMarkerSize(2);

    TH1D* h_unfold3 =(TH1D*)unfold3.Hreco();
    h_unfold3->SetLineWidth(2);
    h_unfold3->SetMarkerSize(2);


    //------Plotting--------------

    //TCanvas* c2 = new TCanvas ("c2", "c2", 1000,1000,1000,800);
    TLegend* leg_hype = new TLegend(0.13,0.7,0.4,0.9);
    Intensity_measured1->SetLineColor(kRed);
    Intensity_measured1->SetLineWidth(1);
    Intensity_measured2->SetLineColor(kBlack);
    Intensity_measured2->SetLineWidth(1);
    Intensity_measured3->SetLineColor(kOrange);
    Intensity_measured3->SetLineWidth(1);
    h->SetLineColor(kBlue);
    h->SetLineWidth(1);
    h_unfold1->SetLineColor(kGreen);
    h_unfold2->SetLineColor(kRed);
    h_unfold3->SetLineColor(kOrange);

    leg_hype->AddEntry(h, "True-data");
    leg_hype->AddEntry(h_unfold1, "Recontruction #sigma = 0.1", "ple");
    leg_hype->AddEntry(h_unfold2, "Recontruction #sigma = 0.3", "ple");
    leg_hype->AddEntry(h_unfold3, "Recontruction #sigma = 0.5", "ple");

    //---------------------------------------------------

    TCanvas* ciccio = new TCanvas("ciccio", "ciccio",1000, 1000,1000, 1000);
    // Upper histogram plot is pad1
    TPad* pad11 = new TPad("pad11", "pad11", 0, 0.28, 1, 1);
    pad11->SetBottomMargin(0);  // joins upper and lower plot
    pad11->SetRightMargin(0.04);
    pad11->SetLeftMargin(0.13);
    pad11->SetTickx(1);
    pad11->SetTicky(1);
    pad11->Draw();
    // Lower ratio plot is pad2
    ciccio->cd();  //returns to main canvas before defining pad2
    TPad* pad22 = new TPad("pad22", "pad22", 0, 0.0, 1, 0.28);
    pad22->SetTopMargin(0); // joins upper and lower plot
    pad22->SetRightMargin(0.04);
    pad22->SetLeftMargin(0.13);
    pad22->SetBottomMargin(0.13);
    pad22->SetGridx();
    pad22->SetGridy();
    pad22->SetTickx(1);
    pad22->SetTicky(1);
    pad22->Draw();

    ciccio->Update();
    pad11->cd();

    gStyle->SetOptStat(0);
    
    h->Draw("hist");
    h->SetTitle("Unfolding result #sigma = {0.1, 0.3, 0.5}");
    h->GetYaxis()->SetTitleSize(.04);
    h->GetXaxis()->SetTitleSize(.04);
    h_unfold1->Draw("SAME");
    h_unfold2->Draw("SAME");
    h_unfold3->Draw("SAME");
    leg_hype->Draw();

    TH1D* rateo1 = createRatio(h, h_unfold1);
    rateo1->SetLineColor(kGreen);
    TH1D* rateo2 = createRatio(h, h_unfold2);
    rateo2->SetLineColor(kRed);
    TH1D* rateo3 = createRatio(h, h_unfold3);
    rateo3->SetLineColor(kOrange);

    pad22->cd();
    rateo1->Draw("hist");
    rateo2->Draw("hist same");
    rateo3->Draw("hist same");

    ciccio->cd();
    ciccio->Update();
    ciccio->Draw();
    ciccio->SaveAs("./p.pdf");

//--------------------------------------

   Intensity_measured1->SetLineColor(kGreen);
   Intensity_measured2->SetLineColor(kRed);
   Intensity_measured3->SetLineColor(kOrange);
   TCanvas* miao = new TCanvas("m","m",1000, 1000,1000, 800);
   //gStyle->SetOptStat(0000);
   TLegend* legg = new TLegend(0.1,0.7,0.4,0.9);

   legg->AddEntry(h, "Clean-data", "ple");
   legg->AddEntry(Intensity_measured1, "#sigma = 0.1");
   legg->AddEntry(Intensity_measured2, "#sigma = 0.3");
   legg->AddEntry(Intensity_measured3, "#sigma = 0.5");

   //Saving the results
   h->Draw("hist");
   Intensity_measured1->Draw("hist same");
   Intensity_measured2->Draw("hist same");
   Intensity_measured3->Draw("hist same");

   legg->Draw("SAME");
   miao->Draw("");

   miao->SaveAs("./smeared_hyp.pdf");

   return 0;

}



 
