#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooStats/HypoTestInverterOriginal.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/HybridCalculatorOriginal.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/HypoTestPlot.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TLatex.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include <iostream>
#include "TF2.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "RooUniform.h"
#include "RooClassFactory.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "RooNLLVar.h"
#include "RooMCStudy.h"
#include "RooChi2Var.h"
#include "RooAddPdf.h"
#include "RooMinimizer.h"
#include <cmath>
#include "string.h"
#include <iostream>

using namespace RooFit ;
using namespace RooStats;

double myPDF(double *x, double *par){ //to define TF2 obj
  return par[3]*(3./(4.*M_PI))*(0.5*(1.-par[0]) + (0.5)*(3.*par[0]-1)*cos(x[0])*cos(x[0]) - par[1]*sin(x[0])*sin(x[0])*cos(2.*x[1])- sqrt(2.)*par[2]*sin(2.*x[0])*cos(x[1]));
}

int main(){
    
    //defining Real variables according to RooFit
    //Angles physical defined in [0-2*pi]
    RooRealVar  t("t","t",0,M_PI) ;
    RooRealVar  p("p","p", 0, 2*M_PI) ;
    
    //Defining Parameters according to RooFit.
    //Initial value is the 3rd argument, 4th and 5th arg
    //are the range of existance of the parameter.
    RooRealVar alpha("alpha", "alpha", 0.65, 0.6, 0.7);

    RooRealVar beta("beta", "beta", 0.06, 0.01, 0.2);
    RooRealVar gamma("gamma", "gamma", -0.18, -1, 0.);
    RooRealVar scale_uni("scale_uni", "scale_uni", 1.);
    RooRealVar for_comp("for_comp", "for_comp", 0.01, 0.001,1);

    RooUniform uni_t("uniformt", "uniformt", t);
    RooUniform uni_p("uniformp", "uniformp", p);
    
    //UNIFORM PDF
    RooAddPdf unif("unif", "uni_t+uni_p", RooArgList(uni_t, uni_p), RooArgList(scale_uni));

    TH2F* h2_pdf = new TH2F("F(#theta, #phi) uniform distribution", "F(#theta, #phi) uniform distribution", 100, 0, M_PI, 100, 0, 2*M_PI);
    unif.fillHistogram(h2_pdf, RooArgList(t,p));
    h2_pdf->GetXaxis()->SetTitle("#theta [rad]");
    h2_pdf->GetYaxis()->SetTitle("#phi [rad]");
    h2_pdf->GetXaxis()->SetTitleOffset(2);
    h2_pdf->GetYaxis()->SetTitleOffset(2);
    h2_pdf->GetZaxis()->SetTitle("F(#theta, #phi)");
    h2_pdf->GetZaxis()->SetTitleOffset(1.5);
    TCanvas* c_pdf = new TCanvas("c_pdf", "c_pdf", 1000,1000,1000,800);
    c_pdf->SetLeftMargin(2);
    h2_pdf->Draw("surf");
    c_pdf->Draw();
    c_pdf->SaveAs("./graphs/real_uniform_pdf.png", "png");
    
    //ClassFactory generates a pdf instance and normalize it.
    RooAbsPdf* genpdf = RooClassFactory::makePdfInstance("GenPdf","(3./(4.*M_PI))*(0.5*(1.-alpha) + (0.5)*(3.*alpha-1)*cos(t)*cos(t) - beta*sin(t)*sin(t)*cos(2.*p)- sqrt(2.)*gamma*sin(2.*t)*cos(p))",RooArgSet(t,p,alpha, beta, gamma)) ;
    RooAddPdf zero_meno("zero_meno", "gen_pdf+uniform", RooArgList(unif,*genpdf), RooArgList(for_comp) );

    zero_meno.Print();

    RooDataSet* data = genpdf->generate(RooArgSet(t,p), 50000);

    RooFitResult* r = zero_meno.fitTo(*data, Save());
    std::cout << "===> Fit Results:" << std::endl;
    r->Print();
    data->Print();



    //--------------------------------------
    //--------Plotting result of fit--------
    //--------------------------------------

    TH2F* dh3 = new TH2F("h3", "h3", 50, 0, M_PI, 50, 0, 2*M_PI);
    data->fillHistogram(dh3, RooArgList(t,p));
    
    TH1* hh_pdf = zero_meno.createHistogram("t,p", 100,100);
    hh_pdf->SetTitle(Form("#splitline{Unbinned ML fit: #alpha %.3f#pm%.3f, #beta %.3f#pm%.3f, #gamma %.3f#pm%.3f,}{for_comp %.3f#pm%.3f}  ", alpha.getVal(),alpha.getError(), beta.getVal(),beta.getError(), gamma.getVal(), gamma.getError(), for_comp.getVal(), for_comp.getError()));
    hh_pdf->GetXaxis()->SetTitle("#theta [rad]");
    hh_pdf->GetXaxis()->SetTitleOffset(2);
    hh_pdf->GetYaxis()->SetTitle("#phi [rad]");
    hh_pdf->GetYaxis()->SetTitleOffset(2);
    hh_pdf->GetZaxis()->SetTitleOffset(2);
    hh_pdf->SetStats(0);
    hh_pdf->SetLineColor(kRed) ;

    TCanvas* c_uni = new TCanvas("c_uni", "c_uni", 1000,1000,1000,800);
    c_uni->SetLeftMargin(0.2);
    gStyle->SetOptStat();
    hh_pdf->Draw("surf3");
    dh3->Draw("same lego");
    c_uni->Draw();
    c_uni->SaveAs("p.png");

    genpdf->Print();
/* 
    //--------------------------------------
    //-----------Statistical Test-----------
    //--------------------------------------

  RooWorkspace*w = new RooWorkspace("w");
  w->import(t);
  w->import(p);
  w->import(zero_meno);
  w->import(*data);

  w->defineSet("obs",RooArgSet(t,p));
  w->defineSet("poi",RooArgSet(for_comp));
  w->defineSet("nuisParams","alpha, beta, gamma, scale");

  ModelConfig zero_model("zero_model", w);
  zero_model.SetPdf(*w->pdf("zero_meno"));
  zero_model.SetObservables(*w->set("obs"));
  zero_model.SetParametersOfInterest(*w->set("poi"));
  w->var("for_comp")->setVal(0.0);  // important only uniform
  zero_model.SetSnapshot(*w->set("poi"));

  // create the alternate (signal+background) ModelConfig with s=50
  ModelConfig uno_model("uno_model", w);
  uno_model.SetPdf(*w->pdf("zero_meno"));
  uno_model.SetObservables(*w->set("obs"));
  uno_model.SetParametersOfInterest(*w->set("poi"));
  w->var("for_comp")->setVal(0.5); // Setting non zero value to genpdf
  uno_model.SetSnapshot(*w->set("poi"));

  //Enforcing pdf on nuisance
  RooAbsPdf * nuisPdf = MakeNuisancePdf(uno_model,"nuisancePdf_uno_model");

  w->Print();


  SimpleLikelihoodRatioTestStat slrts(*zero_model.GetPdf(),*uno_model.GetPdf());
  slrts.SetNullParameters(*zero_model.GetSnapshot());
  slrts.SetAltParameters(*uno_model.GetSnapshot());

  HybridCalculator hc2(*data, zero_model, uno_model);
  ToyMCSampler *toymcs2 = (ToyMCSampler*)hc2.GetTestStatSampler();
  hc2.ForcePriorNuisanceAlt(*nuisPdf);
  hc2.ForcePriorNuisanceNull(*nuisPdf);
  toymcs2->SetTestStatistic(&slrts);
  toymcs2->SetNEventsPerToy(1);
  hc2.SetToys(1000,1000);

  HypoTestResult *r2 = hc2.GetHypoTest();
  std::cout << "-----------------------------------------"<<std::endl;
  std::cout << "Part 5" << std::endl;
  r2->Print();

  //--------------------------------------
  //----------Statistical Test 2----------
  //--------------------------------------


  std::cout<< "--------------------------"<<std::endl;
  std::cout<< "Statistical_Test_2"<<std::endl;
  std::cout<< "--------------------------"<<std::endl;

  FrequentistCalculator *  fc  = new FrequentistCalculator(*data, zero_model, uno_model);
  fc->SetToys(1000,1000);    // 1000 for null (B) and 1000 for alt (S+B) 

  // create the test statistics
  ProfileLikelihoodTestStat profll(*uno_model.GetPdf());
  // use one-sided profile likelihood
  profll.SetOneSidedDiscovery(true);

  // configure  ToyMCSampler and set the test statistics
  ToyMCSampler *toymcs = (ToyMCSampler*)fc->GetTestStatSampler();
  toymcs->SetTestStatistic(&profll);

  if (!uno_model.GetPdf()->canBeExtended())
    toymcs->SetNEventsPerToy(1);

  // run the test
  HypoTestResult * res = fc->GetHypoTest();
  res->Print();

  // plot test statistic distributions
  TCanvas* c_boh = new TCanvas("cboh", "cboh", 1000,1000,1000,800);
  HypoTestPlot * plot = new HypoTestPlot(*res);
  plot->Draw();
  c_boh->Draw();
  c_boh->SaveAs("./boh.png");
  //w.factory("SUM:model(nsig[0,10000]*sig_pdf, nbkg[0,10000]*bkg_pdf)")




    HybridCalculatorOriginal myhc(*data, zero_meno, unif,0,0);
    myhc.SetTestStatistic(2);
    myhc.SetNumberOfToys(1000);
    myhc.UseNuisance(false);

    // run the hypothesis-test invertion
    HypoTestInverterOriginal myInverter(myhc,for_comp);
    
    
    myInverter.SetTestSize(0.10);
    myInverter.UseCLs(true);
     myInverter.RunFixedScan(0,0,1);
    // scan for a 95% UL
    myInverter.RunAutoScan(0,1,myInverter.Size()/2,0.005);
    // run an alternative autoscan algorithm
    //myInverter.RunAutoScan(0,6,myInverter.Size()/2,0.005,1);
    //myInverter.RunOnePoint(3.9);


    HypoTestInverterResult* results = myInverter.GetInterval();

    TCanvas* c2 = new TCanvas("c2", "c2", 1000,1000,1000,800);
    HypoTestInverterPlot myInverterPlot("myInverterPlot","",results);
    TGraphErrors* gr1 = myInverterPlot.MakePlot();
    gr1->Draw("ALP");
    c2->Draw();
    c2->SaveAs("./test.png");

*/

//--------------------------------------
//----------Statistical Test 3----------
//--------------------------------------

/* 

std::cout<< "--------------------------"<<std::endl;
std::cout<< "Statistical_Test_3"<<std::endl;
std::cout<< "--------------------------"<<std::endl;

RooWorkspace*w = new RooWorkspace("w");
  w->import(t);
  w->import(p);
  w->import(*genpdf);
  w->import(*data);

  w->defineSet("obs",RooArgSet(t,p));
  w->defineSet("arg",RooArgSet(alpha, beta, gamma, scale));
  //w->defineSet("nuisParams","alpha, beta, gamma, scale");

  ModelConfig zero_model("zero_model", w);
  zero_model.SetObservables(*w->set("obs"));
  zero_model.SetPdf(*w->pdf("genpdf"));
  

*/




  /* 
  zero_model.SetParametersOfInterest(*w->set("arg"));
  
  w->var("alpha")->setVal(1./3.);  // important only uniform
  w->var("beta")->setVal(0.0);
  w->var("gamma")->setVal(0.0);
  w->var("scale")->setVal(1);
  zero_model.SetSnapshot(*w->set("arg"));

  // create the alternate meaning uno boson.
  ModelConfig uno_model("uno_model", w);
  uno_model.SetPdf(*w->pdf("genpdf"));
  uno_model.SetObservables(*w->set("obs"));
  uno_model.SetParametersOfInterest(*w->set("arg"));
  w->var("alpha")->setVal(0.65);  // important for distribution genpdf
  w->var("beta")->setVal(0.06);
  w->var("gamma")->setVal(-0.18);
  w->var("scale")->setVal(5);
  uno_model.SetSnapshot(*w->set("arg"));



  //Enforcing pdf on nuisance
  RooAbsPdf * nuisPdf = MakeNuisancePdf(uno_model,"nuisancePdf_uno_model");

  w->Print();
 */
    return 0;
}
